#!/usr/bin/env python

'''Calculates Coverage Per Contig and Gene from MagicBlast Tabular Output.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Magic Blast output should be filtered prior to using this script   !!
!! Use MagicBlast01_ShortRead_Filter.py or other method.              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This script calculates Gene coverage based on CDS positional       !!
!! information as well as inter-gene or between CDS regions. This     !!
!! script assumes the CDS file to be ordered by genome position.      !!
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


def calc_genome_coverage(tbf, rgf_tad, rgf_ani, thrshld):
    """ Reads tabblast file and adds coverage by genome position """

    with open(tbf, 'r') as f:
        for l in f:
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

            if pident >= thrshld:
                for i in range(strt, stp+1, 1):
                    rgf_tad[contig_name][i] += 1
                    rgf_ani[contig_name][i].append(pident)

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


def get_strt_stp(location):
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


def retrieve_gene_coverage(pgf, rgf_tad, rgf_ani):
    """ Retrieves list of depths for each bp position of gene length """

    gn_tad = defaultdict(list) # initialize dictionary
    gn_ani = defaultdict(list)
    gn_len = {}

    intergn_tad = defaultdict(list) # initialize dictionary
    intergn_ani = defaultdict(list)
    intergn_len = {}
    intergn_count = 1

    with open(pgf, 'r') as f:
        stp = 1
        for name, seq in read_fasta(f):
            contig_name = '_'.join(name.split('|')[1].split('_')[:2])
            protein = name.split('protein=')[1].split(']')[0]
            try: protein_id = name.split('protein_id=')[1].split(']')[0]
            except: protein_id = 'pseudo-gene'
            locus_tag = name.split('locus_tag=')[1].split(']')[0]
            location = name.split('location=')[1].split(']')[0].split('(')

            p1, p2 = get_strt_stp(location)

            strt = min(p1, p2) # start of CDS region

            # Define intergenic or between CDS regions
            intergene_strt = stp # start of inter-CDS region
            intergene_stp = strt # stop of inter-CDS region
            intergene_len = intergene_stp - intergene_strt + 1
            intergene_name = f'intergene_{intergn_count:06}'
            intergn_count += 1

            stp = max(p1, p2) # stop of CDS region

            gene_name = f'{contig_name}:{locus_tag}:{protein_id}:{protein}'

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
            f'file extension of one of {fqtype} or {fatype}.'
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


def calc_tad_anir_relabndc(
                            mtg,
                            wglen,
                            rgf_tad,
                            rgf_ani,
                            rgf_len,
                            gn_tad,
                            gn_ani,
                            gn_len,
                            intergn_tad,
                            intergn_ani,
                            intergn_len,
                            tad,
                            outpre,
                            precision
                            ):

    """ Calculate tad and anir for whole genome, contig, and gene """

    print('... Calculating TADs for Contigs')
    contig_tad, contig_breadth, wg_tad = get_contig_tad(rgf_tad, tad)

    print('... Calculating TADs for Genes')
    gene_tad, gene_breadth = get_gene_tad(gn_tad, tad)

    print('... Calculating TADs for Inter-Gene Regions')
    intergene_tad, intergene_breadth = get_gene_tad(intergn_tad, tad)

    print('... Calculating TAD for Genome')
    wgbreadth = sum(i > 0 for i in wg_tad) / len(wg_tad)
    wgtad = get_average(wg_tad, tad)

    print('... Calculating Total Metagenome Size & Relative Abundance')
    relabndc, total_metagenome_bp = get_relative_abundance(wg_tad, mtg)

    print('... Calculating ANI for Contigs')
    contig_ani, wg_ani = get_contig_anir(rgf_ani, tad)

    print('... Calculating ANI for Genes')
    gene_ani = get_gene_anir(gn_ani, tad)

    print('... Calculating ANI for Inter-Gene Regions')
    intergene_ani = get_gene_anir(intergn_ani, tad)

    print('... Calculating ANI for Genome')
    wgani = get_average(wg_ani, tad)

    print('\nWriting output files ...')
    _ = write_file(contig_tad, rgf_len, outpre, '_contig_tad.tsv', precision)
    _ = write_file(
        contig_breadth, rgf_len, outpre, '_contig_breadth.tsv', precision
        )
    _ = write_file(contig_ani, rgf_len, outpre, '_contig_anir.tsv', precision)
    _ = write_file(gene_tad, gn_len, outpre, '_gene_tad.tsv', precision)
    _ = write_file(gene_breadth, gn_len, outpre, '_gene_breadth.tsv', precision)
    _ = write_file(gene_ani, gn_len, outpre, '_gene_anir.tsv', precision)
    _ = write_file(
        intergene_tad, intergn_len, outpre, '_inter-gene_tad.tsv', precision
        )
    _ = write_file(
        intergene_breadth, intergn_len, outpre,
        '_inter-gene_breadth.tsv', precision
        )
    _ = write_file(
        intergene_ani, intergn_len, outpre, '_inter-gene_anir.tsv', precision
        )

    return wgtad, wgbreadth, wgani, relabndc, total_metagenome_bp


def operator(mtg, rgf, tbf, pgf, thd, tad, outpre):
    """ Runs the different functions and writes out results """

    tadp = tad / 100
    precision = 2 # number of decimals places to keep.

    print(f'Using values {tad}% for TAD & {thd}% for ANIr')

    print('Preparing base pair array for each contig in genome.')
    rgf_tad, rgf_ani, rgf_len, wglen = read_genome_lengths(rgf)

    print(
        'Calculating coverage for each base pair position in genome.'
        'This can take a while depending on the number of blast results.'
        )
    rgf_tad, rgf_ani = calc_genome_coverage(tbf, rgf_tad, rgf_ani, thd)

    print('Writing Whole genome per base pair depth')
    _ = write_genome_cov_by_bp(rgf_tad, outpre)

    print('Retrieving coverage for each contig & gene')
    (
        gn_tad,
        gn_ani,
        gn_len,
        intergn_tad,
        intergn_ani,
        intergn_len
            ) = retrieve_gene_coverage(pgf, rgf_tad, rgf_ani)

    print(f'Calculating {tad}% truncated average depth and {thd}% ANIr')
    (
        wgtad,
        wgbreadth,
        wgani,
        relabndc,
        total_metagenome_bp
            ) = calc_tad_anir_relabndc(
                                        mtg,
                                        wglen,
                                        rgf_tad,
                                        rgf_ani,
                                        rgf_len,
                                        gn_tad,
                                        gn_ani,
                                        gn_len,
                                        intergn_tad,
                                        intergn_ani,
                                        intergn_len,
                                        tadp,
                                        outpre,
                                        precision
                                        )

    wg_header = (
            f'Genome_Name\tTAD_{int(tad)}\tBreadth\tANIr_{int(thd)}\t'
            f'Relative_Abundance(%)\tGenome_Length(bp)\tMetagenome_Length(bp)\n'
            )
    wg_lineout = (
            f'{outpre}\t{wgtad:.{precision}f}\t{wgbreadth:.{precision}f}\t'
            f'{wgani:.{precision}f}%\t'
            f'{relabndc:.{precision}f}%\t{wglen}\t{total_metagenome_bp}\n'
            )

    with open(f'{outpre}_genome.tsv', 'w') as o:
        o.write(wg_header)
        o.write(wg_lineout)

    print('\nScript seems to have finished successfully.\n')

    print('\nWhole Genome Values:\n')
    print(wg_header[:-1])
    print(wg_lineout[:-1])


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-m', '--metagenome_file',
        help='Please specify the query metagenome fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-g', '--ref_genome_file',
        help='Please specify the genome fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-b', '--tabular_blast_file',
        help='Please specify the tabular blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--CDS_from_genomic_file',
        help='Please specify the prodigal gene fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--pIdent_threshold_cutoff',
        help='Please specify pIdent threshold to use! (ie: 95)',
        metavar='',
        type=float,
        required=True
        )
    parser.add_argument(
        '-d', '--truncated_avg_depth_value',
        help='Please specify TAD value! (ie: 80 or 90)',
        metavar='',
        type=float,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file_prefix',
        help='What do you like the output file prefix to be?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    operator(
            args['metagenome_file'],
            args['ref_genome_file'],
            args['tabular_blast_file'],
            args['CDS_from_genomic_file'],
            args['pIdent_threshold_cutoff'],
            args['truncated_avg_depth_value'],
            args['out_file_prefix']
            )


if __name__ == "__main__":
    main()
