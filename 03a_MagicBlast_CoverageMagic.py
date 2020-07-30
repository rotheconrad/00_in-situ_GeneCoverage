#!/usr/bin/env python

'''Calculates ANIr and Sequence Coverage from Magic Blast Output.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Magic Blast output should be filtered prior to using this script   !!
!! Use 01c_ShortRead_Filter.py or other method.                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This script calculates ANIr and sequence coverage (as depth and    !!
!! breadth) from tabular Magic Blast output for the whole genome or   !!
!! or MAG, each contig in the genomic sequence fasta file, each       !!
!! predicted coding sequence (protein coding gene), and each          !!
!! intergenic region. It requires the metagenomic fasta file used as  !!
!! the blast queries, the genonimic reference sequence used as the    !!
!! database or subject, and a fasta file of predicted protein coding  !!
!! sequences predicted by Prodigal from the genomic reference.        !!
!! The *_CDS_from_genomic.fna files from NCBI are usually in order.   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Coverage calculated as Truncated Average Depth (TAD):
    * Set TAD to 100 for no truncatation.
    * TAD 80 removes the top 10% and bottom 10% of base pair depths and
      caluclates coverage from the middle 80% of values. Intended to 
      reduce effects of conserved motif peaks and contig edge valleys.
    * Coverage = base pairs recruited / length of genome, contig, or gene

Coverage calculated as Breadth:
    * number of positions in reference sequence covered by at least
      one read alignment divided the length of the reference sequence.

Relative Abundance is calculated as:
    * base pairs recruited / base pairs in metagenome * 100
    * It is the percent of base pairs recruited out of the total
      base pairs sequenced in the metagenome.

ANIr is calculated as:
    * average percent identity of sequence alignments for all reads 
      (should be 1 blast match per read)

This tool takes the following input parameters:

    * Tabular Blast file containing results for 1 genome and 1 metagenome
    * *_genomic.fna file from NCBI used as reference for blast search.
    * Metagenome fastq (or fasta) file used as queries for blast search.
    * *_CDS_from_genomic.fna file of predicted genes for *_genomic.fna.

This script returns the following files:

    * 3 column tsv output of Contig(or gene_name), coverage(or ANIr), 
      sequence length.
    * Writes 11 files total:
        - {out_file_prefix}_genome_by_bp.tsv
        - {out_file_prefix}_genome.tsv
        - {out_file_prefix}_contig_tad.tsv
        - {out_file_prefix}_contig_breadth.tsv
        - {out_file_prefix}_contig_anir.tsv
        - {out_file_prefix}_gene_tad.tsv
        - {out_file_prefix}_gene_breadth.tsv
        - {out_file_prefix}_gene_anir.tsv
        - {out_file_prefix}_intergene_tad.tsv
        - {out_file_prefix}_intergene_breadth.tsv
        - {out_file_prefix}_intergene_anir.tsv

*_gene_* files contain values for the CDS regions.
*_intergene_* files contain values for the inter-CDS regions.

This script requires the following packages:

    * argparse, sys
    * collection.defaultdict
    * itertools

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 12th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys
import itertools
from collections import defaultdict


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def read_genome_lengths(rgf):
    """ Reads genome lengths file returns dict genome_name: length """

    rgf_tad = defaultdict(dict) # initialize dicts
    rgf_ani = defaultdict(dict)
    rgf_len = {}
    wg_len = 0

    # read through genome fasta and build dictionary of dictionary
    # Containing base pair position for length of each contig.
    with open(rgf, 'r') as f:
        for name, seq in read_fasta(f):

            contig_name = name.split(' ')[0][1:]
            length = len(seq) # calculate length of contig

            rgf_len[contig_name] = length
            wg_len += length

            # This populates the dictionary with value of zero for each
            # base pair position for each contig in the genome fasta
            for i in range(1, length+1, 1):
                rgf_tad[contig_name][i] = 0
                rgf_ani[contig_name][i] = []

    return rgf_tad, rgf_ani, rgf_len, wg_len


def calc_genome_coverage(tbf, rgf_tad, rgf_ani, lthd, uthd):
    """ Reads tabblast file and adds coverage by genome position """

    read_count = 0
    with open(tbf, 'r') as f:
        for l in f:
            # Progress Tracker
            if read_count % 50000 == 0:
                print(f'... Blast matches processed so far ... {read_count:013}')
            read_count += 1

            # Skip magic blast header
            if l.startswith('#'): continue

            # split each line and define variables of interest
            X = l.rstrip().split('\t')
            pident = float(X[2])
            contig_name = X[1]
            strt = min(int(X[8]), int(X[9]))
            stp = max(int(X[8]), int(X[9]))
            
            # for each read above the user specified threshold, add
            # coverage of +1 to each basepair position for length of
            # read alignment along the subject sequence.

            if pident > lthd and pident < uthd:
                for i in range(strt, stp+1, 1):
                    rgf_tad[contig_name][i] += 1
                    rgf_ani[contig_name][i].append(pident)

    print(f'... Total Blast matches process: {read_count:013}')

    return rgf_tad, rgf_ani


def write_genome_cov_by_bp(rgf_tad, outpre):
    """ writes the tad coverage per bp as 2 col tsv file pos: value """

    counter = 1
    with open(f'{outpre}_genome_by_bp.tsv', 'w') as o:
        o.write('Position\tDepth\n')
        for k, v in rgf_tad.items():
            for i in v.values():
                o.write(f'{counter}\t{i}\n')
                counter += 1


def get_ncbi_strt_stp(location):
    """ Sorts out the NCBI CDS from genomic location nonsense """

    if len(location) == 1:
        X = location[0].split('..')
        p1 = int(''.join(i for i in X[0] if i.isdigit()))
        p2 = int(''.join(i for i in X[1] if i.isdigit()))

    elif len(location) == 2:
        X = location[1].split('..')
        if len(X) == 2:
            p1 = int(''.join(i for i in X[0] if i.isdigit()))
            p2 = int(''.join(i for i in X[1] if i.isdigit()))
        elif len(X) == 3:
            p1 = int(''.join(i for i in X[0] if i.isdigit()))
            p2 = int(''.join(i for i in X[2] if i.isdigit()))
        else:
            print('Gene location error 1 for entry:')
            print(name)
            print(X)
            sys.exit()

    elif len(location) == 3:
        X = location[2].split('..')
        if len(X) == 3:
            p1 = int(''.join(i for i in X[0] if i.isdigit()))
            p2 = int(''.join(i for i in X[2] if i.isdigit()))
        else:
            print('Gene location error 2 for entry:')
            print(name)
            print(X)
            sys.exit()

    else:
        print('Gene location error 3 for entry:')
        print(name)
        print(location)
        sys.exit()

    return p1, p2


def retrieve_ncbi_gene_coverage(pgf, rgf_tad, rgf_ani):
    """ Retrieves list of depths for each bp position of gene length """

    gn_tad = defaultdict(list) # initialize dictionary
    gn_ani = defaultdict(list)
    gn_len = {}

    intergn_tad = defaultdict(list) # initialize dictionary
    intergn_ani = defaultdict(list)
    intergn_len = {}

    with open(pgf, 'r') as f:
        stp = 0
        for name, seq in read_fasta(f):
            contig_name = '_'.join(name.split('|')[1].split('_')[:2])
            try: protein = name.split('protein=')[1].split(']')[0]
            except: protein = 'n/a'
            try: protein_id = name.split('protein_id=')[1].split(']')[0]
            except: protein_id = 'pseudo-gene'
            locus_tag = name.split('locus_tag=')[1].split(']')[0]
            location = name.split('location=')[1].split(']')[0].split('(')

            p1, p2 = get_ncbi_strt_stp(location)

            strt = min(p1, p2) # start of CDS region

            # Define intergenic or between CDS regions
            intergene_strt = stp+1 # start of inter-CDS region
            intergene_stp = strt-1 # stop of inter-CDS region
            intergene_len = intergene_stp - intergene_strt

            stp = max(p1, p2) # stop of CDS region

            gene_name = f'{contig_name}:{locus_tag}:{protein_id}:{protein}'

            intergene_name = (
                f'{contig_name}_intergene_{intergene_strt}-{intergene_stp}'
                )

            gn_len[gene_name] = len(seq)
            intergn_len[intergene_name] = intergene_len

            # Get depth values for gene (CDS) regions
            for i in range(strt, stp+1, 1):
                gn_tad[gene_name].append(rgf_tad[contig_name][i])
                gn_ani[gene_name].extend(rgf_ani[contig_name][i])

            # Get depth values for intergene (inter-CDS) regions
            for i in range(intergene_strt, intergene_stp+1, 1):
                intergn_tad[intergene_name].append(rgf_tad[contig_name][i])
                intergn_ani[intergene_name].extend(rgf_ani[contig_name][i])

        # Get intergene region after last predicted coding region.
        intergene_strt = stp + 1
        intergene_stp = len(rgf_tad[contig_name])
        intergene_len = intergene_stp - intergene_strt

        intergene_name = (
            f'{contig_name}_intergene_{intergene_strt}-{intergene_stp}'
            )
        intergn_len[intergene_name] = intergene_len
        # Get depth values for intergene (inter-CDS) regions
        for i in range(intergene_strt, intergene_stp+1, 1):
            intergn_tad[intergene_name].append(rgf_tad[contig_name][i])
            intergn_ani[intergene_name].extend(rgf_ani[contig_name][i])

    return gn_tad, gn_ani, gn_len, intergn_tad, intergn_ani, intergn_len


def retrieve_prodigal_gene_coverage(pgf, rgf_tad, rgf_ani):
    """ Retrieves list of depths for each bp position of gene length """

    gn_tad = defaultdict(list) # initialize dictionary
    gn_ani = defaultdict(list)
    gn_len = {}

    intergn_tad = defaultdict(list) # initialize dictionary
    intergn_ani = defaultdict(list)
    intergn_len = {}

    with open(pgf, 'r') as f:
        stp = 1 # initial stp value for intergene calculation
        for name, seq in read_fasta(f):
            X = name.split(' # ')
            gene_name = X[0][1:]
            contig_name = '_'.join(gene_name.split('_')[:-1])

            strt = min(int(X[1]), int(X[2]))

            # Define intergenic or between CDS regions
            intergene_strt = stp + 1 # start of inter-CDS region
            intergene_stp = strt # stop of inter-CDS region
            intergene_len = intergene_stp - intergene_strt

            stp = max(int(X[1]), int(X[2]))

            intergene_name = (
                f'{contig_name}_intergene_{intergene_strt}-{intergene_stp}'
                )

            gn_len[gene_name] = len(seq)
            intergn_len[intergene_name] = intergene_len

            # Get depth values for gene (CDS) regions
            for i in range(strt, stp+1, 1):
                gn_tad[gene_name].append(rgf_tad[contig_name][i])
                gn_ani[gene_name].extend(rgf_ani[contig_name][i])

            # Get depth values for intergene (inter-CDS) regions
            for i in range(intergene_strt, intergene_stp+1, 1):
                intergn_tad[intergene_name].append(rgf_tad[contig_name][i])
                intergn_ani[intergene_name].extend(rgf_ani[contig_name][i])

        # Get intergene region after last predicted coding region.
        intergene_strt = stp + 1
        intergene_stp = len(rgf_tad[contig_name])
        intergene_len = intergene_stp - intergene_strt

        intergene_name = (
            f'{contig_name}_intergene_{intergene_strt}-{intergene_stp}'
            )
        intergn_len[intergene_name] = intergene_len
        # Get depth values for intergene (inter-CDS) regions
        for i in range(intergene_strt, intergene_stp+1, 1):
            intergn_tad[intergene_name].append(rgf_tad[contig_name][i])
            intergn_ani[intergene_name].extend(rgf_ani[contig_name][i])

    return gn_tad, gn_ani, gn_len, intergn_tad, intergn_ani, intergn_len


def truncate(x, tad):
    """ returns tad range of a list/array """

    xsorted = sorted(x)
    xlen = len(x) # get length of the list
    inverse_tad = round((1.0 - tad) / 2.0, 2) # to get top and bottom 
    q = int(xlen * inverse_tad) # to get top and bottom
    bottom = q
    top = xlen - q
    tad_range = xsorted[bottom:top] # slice list

    #print(xlen, bottom, top, len(tad_range), sum(xsorted[:bottom+1]))

    return tad_range


def get_contig_tad(rgf_tad, tad):
    """ reads through rgf_tad and returns dict of tads by contig """

    contig_tad = {}
    contig_breadth = {}
    wg_tad = []

    for k, v in rgf_tad.items():
        values = list(v.values())
        breadth = sum(i > 0 for i in values) / len(values)
        contig_breadth[k] = breadth
        coverage = get_average(values, tad)
        contig_tad[k] = coverage
        wg_tad.extend(values)

    return contig_tad, contig_breadth, wg_tad


def get_gene_tad(gn_tad, tad):
    """ reads through gn_tad and returns dict of tads by gene """

    gene_tad = {}
    gene_breadth = {}

    for k, v in gn_tad.items():
        breadth = sum(i > 0 for i in v) / len(v)
        gene_breadth[k] = breadth
        gene_tad[k] = get_average(v, tad)

    return gene_tad, gene_breadth
    

def get_contig_anir(rgf_ani, tad):
    """ loops through in_d and calculates tad/ani for each key """

    contig_ani = {}
    wg_ani = []

    for k, v in rgf_ani.items():
        values = list(itertools.chain.from_iterable(list(v.values())))
        average = get_average(values, tad)
        if average > 0: contig_ani[k] = average
        wg_ani.extend(values)

    return contig_ani, wg_ani


def get_gene_anir(gn_ani, tad):
    """ loops through in_d and calculates tad/ani for each key """

    gene_ani = {}

    for k, v in gn_ani.items():
        average = get_average(v, tad)
        if average > 0: gene_ani[k] = average

    return gene_ani


def get_average(l, tad):
    """ returns anir from truncated list """

    trunc_val = truncate(l, tad)

    if sum(trunc_val) == 0:
        average = 0
    else:
        average = sum(trunc_val) / len(trunc_val)

    return average


def get_relative_abundance(wg_tad, mtg):
    """ calculates and returns relative abundance from wg_TAD """

    total_metagenome_bp = 0

    # check if fasta or fastq
    file_type = mtg.split('.')[-1]
    fqtype = ['fastq', 'fq']
    fatype = ['fasta', 'fna', 'fst', 'fa']

    if file_type in fqtype:
        line_count = 0

        with open(mtg, 'r') as f:
            for l in f:
                line_count += 1
                if line_count%4 == 0:
                    total_metagenome_bp += len(l.rstrip())
    elif file_type in fatype:
        with open(mtg, 'r') as f:
            for name, seq in read_fasta(f):
                total_metagenome_bp += len(seq)
    else:
        print(
            'Error in determining metagenome format of fasta or fastq. '
            'Please double check metagenome file type and try again. '
            'Metagenome file should be either fasta or fastq format with '
            f'file extension of one of {fqtype} or {fatype}.\n\n'
            'If there is a file extension you would like added, please '
            'submit a feature request through the issues tab of the '
            'GitHub repo at: '
            'https://github.com/rotheconrad/00_in-situ_GeneCoverage/issues'
            '.\n I will be happy to add aditional file extensions.\n\n'
            )
        sys.exit()

    relabndc = (sum(wg_tad) / total_metagenome_bp) * 100

    return relabndc, total_metagenome_bp


def write_file(in_d, len_d, outpre, outpost, precision):
    """ writes dictionary to file """
    
    outfile = outpre + outpost

    with open(outfile, 'w') as o:
        o.write('Name\tValue\tLength\n')
        for k, v in in_d.items():
            o.write(f'{k}\t{v:.{precision}f}\t{len_d[k]}\n')


def calc_contig_stats(
    rgf_tad, rgf_ani, rgf_len, tad, outpre, precision
    ):
    """Calculate ANIr, TAD and breadth and write to files for Contigs"""

    print('... Calculating TADs for Contigs.')
    contig_tad, contig_breadth, wg_tad = get_contig_tad(rgf_tad, tad)

    print('... Writing Contig TAD file.')
    _ = write_file(contig_tad, rgf_len, outpre, '_contig_tad.tsv', precision)

    print('... Writing Contig Breadth file.')
    _ = write_file(
        contig_breadth, rgf_len, outpre, '_contig_breadth.tsv', precision
        )

    contig_tad = None
    contig_breadth = None

    print('... Calculating ANIr for Contigs')
    contig_anir, wg_ani = get_contig_anir(rgf_ani, tad)

    print('... Writing Contig ANIr file.')
    _ = write_file(contig_anir, rgf_len, outpre, '_contig_anir.tsv', precision)

    contig_anir = None

    return wg_tad, wg_ani


def calc_gene_stats(
    gn_tad, gn_anir, gn_len, tad, outpre, precision
    ):
    """Calculate ANIr, TAD and breadth and write to files for Genes"""

    print('... Calculating TADs for Genes')
    gene_tad, gene_breadth = get_gene_tad(gn_tad, tad)

    print('... Writing Gene TAD file.')
    _ = write_file(gene_tad, gn_len, outpre, '_gene_tad.tsv', precision)

    print('... Writing Gene Breadth file.')
    _ = write_file(gene_breadth, gn_len, outpre, '_gene_breadth.tsv', precision)

    gene_tad = None
    gene_breadth = None

    print('... Calculating ANI for Genes')
    gene_anir = get_gene_anir(gn_anir, tad)

    print('... Writing Gene ANIr file.')
    _ = write_file(gene_anir, gn_len, outpre, '_gene_anir.tsv', precision)

    gene_anir = None


def calc_intergene_stats(
    intergn_tad, intergn_anir, intergn_len, tad, outpre, precision
    ):

    """Calculate ANIr, TAD and breadth and write to Intergene files"""

    print('... Calculating TADs for Inter-Gene Regions')
    intergene_tad, intergene_breadth = get_gene_tad(intergn_tad, tad)

    print('... Writing Intergene TAD file.')
    _ = write_file(
        intergene_tad, intergn_len, outpre, '_inter-gene_tad.tsv', precision
        )

    print('... Writing Intergene Breadth file.')
    _ = write_file(
        intergene_breadth, intergn_len, outpre,
        '_inter-gene_breadth.tsv', precision
        )

    intergene_tad = None
    intergene_breadth = None

    print('... Calculating ANI for Inter-Gene Regions')
    intergene_anir = get_gene_anir(intergn_anir, tad)

    print('... Writing Intergene ANIr file.')
    _ = write_file(
        intergene_anir, intergn_len, outpre, '_inter-gene_anir.tsv', precision
        )

    intergene_anir = None


def calc_genome_stats(
    mtg, wg_tad, wg_anir, wglen, tad, lthd, uthd, outpre, precision
    ):
    """Calculate ANIr, TAD and breadth and writes to Genome files"""

    print('... Calculating TAD for Genome')
    wgbreadth = sum(i > 0 for i in wg_tad) / len(wg_tad)
    wgtad = get_average(wg_tad, tad)

    if mtg is None:
        relabndc = 'n/a'
        total_metagenome_bp = 'n/a'
    elif mtg.isdigit():
        total_metagenome_bp = int(mtg)
        relabndc = (sum(wg_tad) / total_metagenome_bp) * 100
    else:
        print('... Calculating Total Metagenome Size & Relative Abundance')
        relabndc, total_metagenome_bp = get_relative_abundance(wg_tad, mtg)
        relabndc = f'{relabndc:.{precision}f}%'

    print('... Calculating ANI for Genome')
    wganir = get_average(wg_anir, tad)

    wg_header = (
            f'Genome_Name\tTAD_{tad*100}\tBreadth\tANIr_{lthd}-{uthd}\t'
            f'Relative_Abundance(%)\tGenome_Length(bp)\tMetagenome_Length(bp)\n'
            )
    wg_lineout = (
            f'{outpre}\t{wgtad:.{precision}f}\t{wgbreadth:.{precision}f}\t'
            f'{wganir:.{precision}f}%\t'
            f'{relabndc}\t{wglen}\t{total_metagenome_bp}\n'
            )

    with open(f'{outpre}_genome.tsv', 'w') as o:
        o.write(wg_header)
        o.write(wg_lineout)

    print('\nWhole Genome Values:\n')
    print(wg_header[:-1])
    print(wg_lineout[:-1], '\n')


def operator(
    mtg, rgf, tbf, pgf, lthd, uthd, tad, outpre, ncbi
    ):
    """ Runs the different functions and writes out results """

    tadp = tad / 100
    precision = 2 # number of decimals places to keep.

    print(
        f'Using values {tad}% for TAD & {lthd}%, {uthd}% for pIdent Threshold'
        )

    print('\nPreparing base pair array for each contig in genome.')
    rgf_tad, rgf_anir, rgf_len, wglen = read_genome_lengths(rgf)

    print(
        '\nCalculating coverage for each base pair position in genome. \n'
        'This can take a while depending on the number of blast results.'
        )
    rgf_tad, rgf_anir = calc_genome_coverage(tbf, rgf_tad, rgf_anir, lthd, uthd)

    print('\nWriting read depth per base pair to file.')
    _ = write_genome_cov_by_bp(rgf_tad, outpre)

    ### Check for Prodigal or NCBI. #####################################
    if pgf:
        print(
            '\nUsing Prodigal protein fasta.\n'
            'Retrieving coverage for each contig & gene.'
            )

        (
        gn_tad, gn_anir, gn_len, intergn_tad, intergn_anir, intergn_len
            ) = retrieve_prodigal_gene_coverage(pgf, rgf_tad, rgf_anir)

        do_genes = True

    elif ncbi:
        print(
            '\nUsing NCBI CDS from genomic FASTA.\n'
            'Retrieving coverage for each contig & gene'
            )

        (
        gn_tad, gn_anir, gn_len, intergn_tad, intergn_anir, intergn_len
            ) = retrieve_ncbi_gene_coverage(ncbi, rgf_tad, rgf_anir)

        do_genes = True

    else:
        print(
            '\n\n!! No gene prediction file entered.\n'
            '!! NOT Calculating values by genes.\n'
            '!! Calculating values for by contig and whole genome only.\n'
            )

        do_genes = False

    ### End Predicted Gene File Check ###################################

    print(
        f'\nCalculating {tad}% truncated average depth '
        f'and {lthd}%, {uthd}% ANIr'
        )

    wg_tad, wg_anir = calc_contig_stats(
                    rgf_tad, rgf_anir, rgf_len, tadp, outpre, precision
                    )

    # Clear from memory
    rgf_tad = None
    rgf_anir = None
    rgf_len = None

    if do_genes == True:

        _ = calc_gene_stats(
            gn_tad, gn_anir, gn_len, tadp, outpre, precision
            )

        # Clear from memory
        gn_tad = None
        gn_anir = None
        gn_len = None

        _ = calc_intergene_stats(
            intergn_tad, intergn_anir, intergn_len, tadp, outpre, precision
            )

        intergn_tad = None
        intergn_anir = None
        intergn_len = None

    _ = calc_genome_stats(
        mtg, wg_tad, wg_anir, wglen, tadp, lthd, uthd, outpre, precision
        )

    print('\nScript seems to have finished successfully.\n')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-m', '--metagenome_file',
        help=
            '(Optional) To calculate relative abundance specify either the path'
            ' to the metagenome file, or the size of the metagenome in base '
            'pairs.',
        #metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-g', '--ref_genome_file',
        help='Please specify the genome fasta file!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-b', '--tabular_blast_file',
        help='Please specify the tabular blast file!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--prodigal_protein_fasta',
        help=
            '(Optional) Use this option to report values for genes and '
            'intergenic regions predicted with Prodigal.',
        #metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-n', '--NCBI_CDS_genomic',
        help=
            '(Optional) Use this option to report values for genes and '
            ' intergenic regions from an NCBI CDS from genomic FASTA file.',
        #metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-c', '--pIdent_threshold_lower',
        help=
            '(Optional) Lower percent sequence identity of reads to include '
            'coverage calculations (Default = 94.99). [pIdent > value].',
        #metavar='',
        type=float,
        required=False,
        default=94.99
        )
    parser.add_argument(
        '-u', '--pIdent_threshold_upper',
        help=
            '(Optional) Upper percent sequence identity of reads to include '
            'coverage calculations (Default = 100.01). [pIdent < value].',
        #metavar='',
        type=float,
        required=False,
        default=100.01
        )
    parser.add_argument(
        '-d', '--truncated_avg_depth_value',
        help='(Optional) Specify a different TAD value! (Default = 80)',
        #metavar='',
        type=float,
        required=False,
        default=80
        )
    parser.add_argument(
        '-o', '--out_file_prefix',
        help='What do you like the output file prefix to be?',
        #metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')
    operator(
            args['metagenome_file'],
            args['ref_genome_file'],
            args['tabular_blast_file'],
            args['prodigal_protein_fasta'],
            args['pIdent_threshold_lower'],
            args['pIdent_threshold_upper'],
            args['truncated_avg_depth_value'],
            args['out_file_prefix'],
            args['NCBI_CDS_genomic']
            )


if __name__ == "__main__":
    main()
