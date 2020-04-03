# Workflow to calculate ANIr and sequence coverage (depth and breadth) of genome(s) / MAG(s) from metagenomes by gene, intergenic region, contig, and whole genome.

This workflow produces separate files in tab separated value (tsv) format for ANIr, sequence depth, and sequence breadth for the genes, intergenic regions, and contigs of a genome / MAG in fasta format. It also produces a file containing sequence depth at each position of the genome as well as a file with results calculated for the whole genome sequence. tsv files can be easily opened in Excel, imported into Python with Pandas, or read into R for further analysis. There is also an option to generate some summary plots or construct a Data Table with whole genome results for multiple genomes / MAGs across many metagenome samples.

*Additionally, this workflow can be used with Genomic FASTA and CDS from genomic FASTA files retrieved from the [NCBI assembly database](https://www.ncbi.nlm.nih.gov/assembly/). In this case, skip the renaming step for sequence names in the reference fasta file in Step 01 and skip all of Step 02. Use the -n flag for NCBI in Step 03.*

All of the Python scripts in this repository are written for Python version 3.6+. They can be executed and help output obtained by entering:

```bash
python scriptname.py -h
```

#### Motivation:

Read mapping is commonly used to calculate sequencing depth or relative abundance for genomic reference sequences of interest. Some tools report a percent sequence identity metric, while other tools do not. Short read aligners were designed for single genome read mapping, not multi-genome, multi-population metagenomic data from natural environments. In these environments, read mapping tools report reads that stretch outside of the sequence-discrete population range and should be considered as false positives. These false positive matches are reads coming from conserved regions or other closely related populations and should not be included when calculating metrics such as sequencing depth. While a straight line is not the ideal cutoff, there is currently no way a priori to determine a complex curve that would most closely resemble the truth.

This protocol uses MagicBlast for read mapping, Prodigal or NCBI CDS from genomic fasta files for gene prediction, and a series of Python scripts to filter off-target reads and calculate coverage metrics. MagicBlast is a fast short read aligner that has an option to include the Blast style percent sequence identity of read alignments in a tabular output file. The percent sequence identity is used to locate the sequence-discrete population threshold value and then the coverage metrics are calculated using only read alignments above or equal to this value. The threshold value can also be adjusted up or down to explore how the coverage calculations change as the cutoff is changed. This threshold may be higher or lower depending on the characteristics of the population. Some heterogenous populations have a wider distribution while some more clonal populations have a much more narrow distribution.

#### Coverage calculated as Truncated Average Depth (TAD):
- TAD 80 removes the top 10% and bottom 10% of base pair depths and caluclates coverage from the middle 80% of values. Intended to reduce effects of conserved motif peaks and contig edge valleys.
- Coverage = base pairs recruited / length of genome, contig, intergenic region, or gene
- The TAD value for the \*_genome.tsv file (which is also printed to the screen while running the script) is calculated by base pair for all positions in the reference genome or MAG. The genes are ignored at this level and the calculation includes inter-genic regions. The sequencing depth of each base pair position is calculated, the values are numerically sorted, and for TAD 80 the top 10% and bottom 10% of sequence depth values are removed. Then, the average sequencing depth is calculated from the remaining base pair positions.
- For the \*_contig_tad.tsv file this same process is repeated but for each contig in the reference genome or MAG instead of the entire genome.
- For the \*_gene_tad.tsv file this same process is repeated but for each genic region in the Prodigal or CDS from genomic FASTA file.
- For the \*_intergene_tad.tsv file this same process is repeated but for each region between sequences outlined in the Prodigal or CDS from genomic FASTA file.
- Set TAD to 100 for no truncatation.

#### Coverage calculated as Breadth:
- The number of positions in the reference sequence covered by at least one read alignment divided by the length of the reference sequence.

#### Relative Abundance is calculated as:
- base pairs recruited / base pairs in metagenome * 100
- It is the percent of base pairs recruited out of the total base pairs sequenced in the metagenome.

#### ANIr is calculated as:
- The average percent identity of sequence alignments for all reads that map above the user specified sequence identity threshold.

#### This workflow leads to the following result files:

- 3 column tsv output of (Contig, Gene, or Intergenic region), (Depth, Breadth, or ANIr), sequence length.
- Writes 11 files total:
    - \{out_file_prefix\}_genome_by_bp.tsv (make optional)
    - \{out_file_prefix\}_genome.tsv
    - \{out_file_prefix\}_contig_tad.tsv 
    - \{out_file_prefix\}_contig_breadth.tsv
    - \{out_file_prefix\}_contig_anir.tsv (make optional)
    - \{out_file_prefix\}_gene_tad.tsv (optional)
    - \{out_file_prefix\}_gene_breadth.tsv (optional)
    - \{out_file_prefix\}_gene_anir.tsv (optional - add additional option)
    - \{out_file_prefix\}_intergene_tad.tsv (optional)
    - \{out_file_prefix\}_intergene_breadth.tsv (optional)
    - \{out_file_prefix\}_intergene_anir.tsv (optional - add additional option)


## Step 00: Required tools :: Python 3.6+, Prodigal and Magic Blast.


### Python 3.6+ for running the Python scripts in this repo.

Information for installing and running Python can be found [here](https://www.python.org/). I recommend installing [mini conda](https://docs.conda.io/en/latest/miniconda.html) first and then creating an environment for Python 3.6+ and other tools for the project at hand.

*All Python scripts in this repo were written for Python 3.6+. If you get a syntax error the first time you run a script, please first check your Python version.*

### Prodigal for protein coding gene prediction.
 
Information and installation instructions for Prodigal can be found [here](https://github.com/hyattpd/Prodigal). The publication is [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/).

Prodigal can also be installed with a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```bash
conda create -n prodigal
conda activate prodigal
conda install -c bioconda prodigal
```

*If you're using files from the NCBI Assembly database you do not need Prodigal.*

### Magic Blast for short read metagenome mapping to reference genome(s) / MAG(s).

Information and installation instructions for Magic Blast can be found [here](https://ncbi.github.io/magicblast/). The publication is [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2996-x). NOTE: If the latest version of Magic Blast gives errors on your system, navigate to the parent directory and try a previous version.


## Step 01: Map metagenomic reads to reference genome(s) / MAG(s).

### Check metagenome read names and rename if needed. (fastq or fasta).

*Magic Blast cuts the query name at the first white space character and reports this as the query ID.*

>Fastq files can be named as:  
>    @D00468:261:HYTMHBCX2:1:1101:9119:31637 1:N:0:CAGAGAGG+ACTGCATA  
>    @D00468:261:HYTMHBCX2:1:1101:9119:31637 2:N:0:CAGAGAGG+ACTGCATA  
>where the unique identifier comes after the white space chacter.

*Filtering for best hit and retrieving the sequence becomes much more complicated and computationally intensive once the unique identifier is lost. Rename fastq or fasta files before running Magic Blast. I typically assign a short unique sample name to my metagenome files (metagenomeID). I then rename the reads using this short metagenomeID and a read number like so: metagenomeID_readNumber. This can be accomplished using either script below depending if your reads are in fastq or fasta format.*

For fastq formatted metagenome read files:

```bash
# To Display the program description and parameter options
python 01a_Fastq_rename_sequences.py -h

# Example execution:
python 01a_Fastq_rename_sequences.py -i metagenome_file_name.fastq -p metagenomeID
```

For fasta formatted metagenome read files:

```bash
# To Display the program description and parameter options
python 01b_Fasta_rename_sequences.py -h

# Example execution:
python 01b_Fasta_rename_sequences.py -i metagenome_file_name.fastq -p metagenomeID
```

### Check sequence names in reference fasta files and rename if needed.

*Magic Blast truncates sequences names at 50 characters. Prodigal appends the predicted gene number to the end of the sequence names of the contigs. Depending on how many genes you have (2000-8000 typical for microbial genome) your sequence names need to leave enough room to retain the gene number. I typically assign short unique genome identifiers (uniqueID) to the file names of my genomic fasta files or MAGs. I then rename the contigs in each fasta file using this short uniqueID and the contig number like so: uniqueID_contigNumber. This can be accomplished with the following script:*

```bash
# To Display the program description and parameter options
python 01b_Fasta_rename_sequences.py -h

# Example execution:
python 01b_Fasta_rename_sequences.py -i genomic_fasta.fna -p uniqueID
```

*genomic_fasta.fna can be complete or draft genomes or a MAG. Any genomic sequence you are using as your reference sequence.*

### For individual genome or MAG:

1. Make Magic Blast database.

    *If you also have Blast+ installed on your system, make certain you are calling the makeblastdb program that comes with Magic Blast and not the version that comes with Blast+. Try: which makeblastdb*

    ```bash
    makeblastdb -dbtype nucl -in Ref_Genome.fasta -out Ref_Genome.fasta -parse_seqids
    ```

    *If you forget the -parse_seqids flag it will cause errors later.*

2. Run Magic Blast.

    *For the outfile_name of the -out flag use the naming scheme of uniqueID_metagenomeID.blast where uniqueID is the unique identifier for your genome or MAG.*

    ```bash
    magicblast -query {metagenome_fasta} -db Ref_Genome.fasta -infmt (fasta or fastq) -no_unaligned -splice F -outfmt tabular -parse_deflines T -out {outfile_name}.blast
    ```

    *If you forget to set the -parse_deflines flag to true it will cause errors later.*

3. ~~Shuffle blast results.~~

    ~~*The Magic Blast results are output in an ordered format. The filter script keeps the first best match which will bias the best match selection of tied results to the first match. Using the blast command shuf will randomize the order of the Magic Blast results file to prevent this bias.*~~

    ```bash
    shuf {outfile_name}.blast > {outfile_name}.shuf.blast
    ```

    **The curent version of the filter script 01c_MagicBlast_ShortRead_Filter.py takes care of the randomization and the shuf step is no longer neccessary if using it. If you are using a different script for filtering you may still need to shuffle the result.**

4. Filter results for best hits.

    *Magic Blast will report multiple results per metagenomic read. For this analysis we only want to count each read once. Magic Blast will also report short sequence alignments of high identity. If a sequence alignment is 20 base pairs but the read is 150 base pairs this is considered to be a wrong match so we remove it. The -pml flag uses a ratio of alignment length / read length to identify results of this type. A value of 0.7, 0.8 or 0.9 is recommended. Defaults are set to 70 bp minium read length and 0.7 alignment length / read length. Use -h to see the help file for parameter usage to change the defaults.*

    ```bash
    # To Display the program description and parameter options
    python 01c_MagicBlast_ShortRead_Filter.py -h

    # Example execution:
    python 01c_MagicBlast_ShortRead_Filter.py -i {outfile_name}.shuf.blast
    ```

### Competitive read recruitment for multiple genomes or MAGs:

1. Append a unique identifier to the beginning of the sequence name for all contigs in the genome or MAG files if you have not done so already.

    Adjst the cut parameters to select a unique ID from your genomic fasta files to append to the begginging of your sequence names for each genome / MAG.

    ```bash
    for file in *.fna
        do
            uniqueID=`basename $file | cut -d _ -f 1
            sed -i "s/>/>${uniqueID}_/g" $file
        done
    ```

2. Concatenate all genomes / MAGs into a single fasta file.

    ```bash
    cat *.fna >> Combined_Genomes.fasta
    ```
    *Change .fna to match the file extention of your genomes / MAGs*

3. Make Magic Blast database.

    *If you also have Blast+ installed on your system, make certain you are calling the makeblastdb program that comes with Magic Blast and not the version that comes with Blast+ Try: which makeblastdb*

    ```bash
    makeblastdb -dbtype nucl -in Combined_Genomes.fasta -out Combined_Genoems.fasta -parse_seqids
    ```

    *If you forget the -parse_seqids flag it will cause errors later*

4. Run Magic Blast.

    Magic Blast has the ability to keep track of [paired reads](https://ncbi.github.io/magicblast/cook/paired.html) using the option -paired or -query_mate which it reports in column 23 of the [tabular output](https://ncbi.github.io/magicblast/doc/output.html). It is possible for a read pair to align to separate genomes or MAGs in a competitive read recruitment which can cause an issue when using grep to de-concatenate like the example in (7) below. However, for coverage calculations it is not necessary to keep track of paired reads. The forward and reverse reads can be combined into a single file and run through Magic Blast as a single -query file. If you run Magic Blast using a paired option, you will need to use awk or some other method to de-concatenate your competitive blast selecting only the 2nd column of the tabular output otherwise information in column 23 can cause errors downstream.

    ```bash
    cat readpair1.fastq readpair2.fast[aq] > combined_reads.fast[aq]
    ```

    *For the outfile_name of the -out flag use the naming scheme of uniqueID_metagenomeID.blast where uniqueID is the unique identifier for your genome or MAG.*

    ```bash
    magicblast -query combined_reads.fast[aq] -db Combined_Genomes.fasta -infmt (fasta or fastq) -no_unaligned -splice F -outfmt tabular -parse_deflines T -out uniqueID_metagenomeID.blast
    ```

    *If you forget to set the -parse_deflines flag to true it will cause errors later.*

5. ~~Shuffle blast results.~~

    ~~*The Magic Blast results are output in an ordered format. The filter script keeps the first best match which will bias the best match selection of tied results to the first match. Using the blast command shuf will randomize the order of the Magic Blast results file to prevent this bias.*~~

    ```bash
    shuf {outfile_name}.blast > {outfile_name}.shuf.blast
    ```

    **The curent version of the filter script 01c_MagicBlast_ShortRead_Filter.py takes care of the randomization and the shuf step is no longer neccessary if using it. If you are using a different script for filtering you may still need to shuffle the result.**

6. Filter results for best hits.

    *Magic Blast will report multiple results per metagenomic read. For this analysis we only want to count each read once. Magic Blast will also report short sequence alignments of high identity. If a sequence alignment is 20 base pairs but the read is 150 base pairs this is considered to be a wrong match so we remove it. The -pml flag uses a ratio of alignment length / read length to identify results of this type. A value of 0.7, 0.8 or 0.9 is recommended. Defaults are set to 70 bp minium read length and 0.7 alignment length / read length. Use -h to see the help file for parameter usage to change the defaults.*

    ```bash
    # To Display the program description and parameter options
    python 01_MagicBlast_ShortRead_Filter.py -h

    # Example execution:
    python 01_MagicBlast_ShortRead_Filter.py -i uniqueID_metagenomeID.shuf.blast
    ```

    *The optional -rtm flag can be useful for calculating unique coverage distances between closely related genomes or MAGs. It will remove all reads with a tied best-hit blast match reducing the coverage calculations for shared (or core) sequence regions and emphasizing the unique sequence components present in-situ.*

7. De-concatenate.

    *Now we want to retrieve the results for each genome or MAG from the concatenated results.*

    First, compile a list of the unique genome or MAG IDs:

    ```bash
    for file in *.fna
        do
            uniqueID=`basename $file | cut -d _ -f 1`
            echo ${uniqueID} >> uniqueID_list.txt
        done
    ```

    Then we iterate through the list of unique IDs and place all the blast matches with that ID in new files:

    ```bash
    while read uniqueID
        do
            grep ${uniqueID} uniqueID_metagenomeID.fltrdBstHts.blst >> ${uniqueID}.blast
        done < uniqueID_list.txt
    ```

    *Another method that works even if you ran Magic Blast with paired option thanks to Nastassia Patin @microbesatsea*

    ```bash
    while read uniqueID
        do
            awk -F '\t' '$2 ~ /$uniqueID/' uniqueID_metagenomeID.fltrdBstHts.blst >> ${uniqueID}.blast
        done < uniqueID_list.txt
    ```

    *And additional options to de-concatenate from Nastassia Patin @microbesatsea*

    ```bash
    #To split one filtered .blst file by a list of MAGs in the file MAGs_list.txt:
    while read line; do
        awk -v line="$line" -F '\t' '$2 ~ line' file.blst >> 'file_'$line'.tsv'
    done < MAGs_list.txt

    #To split several filtered .blst files by the same list, retaining blast file name and MAG name in the output file names:
    for f in *.blst; do
        n=`echo $f | cut -d _ -f 1` &&
            while read line; do 
                awk -v line="$line" -F '\t' '$2 ~ line' $f >> $n'_'$line'.tsv';
            done < MAGs_list.txt
    done
```


## Step 02: Predict protein coding genes with Prodigal.

As per the Prodigal documentation, Prodigal can be easily run like this:

```bash
prodigal -i genomic_fasta.fna -o my.genes -a my.proteins.faa
```

*We only care about the my.proteins.faa file for the purpose of this pipeline but it can also be handy to use the -f and -d flags at the same time to generate additional file types for use later.*


## Step 03: Calculate ANIr and Coverage.

*The script takes 1 uniqueID.blast file at a time with its corresponding metagenome, genomic fasta, and predicted genes files. For the outfile_prefix of the -o flag use the naming scheme of uniqueID_metagenomeID where uniqueID is the unique identifier for your genome or MAG. The default value for the TAD calculations is set to 80. The defualt value for minimum and maximum percent identity of read alignments to use in the calculations are 94.99 and 100.1. Use -h to see the help file for parameter usage to change the defaults.*

If using NCBI Assembly Files:
```bash
# To Display the program description and parameter options
python 03a_MagicBlast_CoverageMagic.py -h

# Example execution:
python 03a_MagicBlast_CoverageMagic.py -m metagenomeID.fna -g uniqueID_genomic_FASTA.fna -n uniqueID_CDS_from_genomic_FASTA.fna -b uniqueID.blast -c 95 -d 80 -o uniqueID_metagenomeID
```

If using Prodigal:
```bash
# To Display the program description and parameter options
python 03a_MagicBlast_CoverageMagic.py -h

# Example execution:
python 03a_MagicBlast_CoverageMagic.py -m metagenomeID.fna -g uniqueID.fna -p my.proteins.faa -b uniqueID.blast -c 95 -d 80 -o uniqueID_metagenomeID
```

*If using Genomic FASTA and CDS from genomic FASTA files retrieved from the NCBI assembly database, the Genomic FASTA goes to the -g flag and CDS from genomic FASTA goes to the -n flag. If using Prodigal, the my.proteins.faa file goes to the -p flag and the genomic reference goes to the -g flag.*

The -d flag is for the truncated average value or TAD parameter. A value of 100 will return results without truncation.

The -c flag is a cutoff threshold for the percent identity of the metagenomic read alignments to the genomic reference. Sequence discontinuity gaps for sequence discrete populations are generally observed around 95% percent sequence identity. This is a good starting point, but depending on your target population you may want to increase or decrease this value. Looking at the distribution of percent identity values or a recruitment plot is a great way to investigate the sequence discontinuity for your population of interest.

You can visualize your sequence identity distribution with a histogram using the following:

```bash
# To Display the program description and parameter options
python 03b_MagicBlast_pIdent_Hist.py -h

# Example execution:
python 03b_MagicBlast_pIdent_Hist.py -i uniqueID.blast
```

Example Histogram Plot:

![alt text](03c_Example_plot.png "Example histogram plot.")

Or you can use the [Enveomics collection](http://enve-omics.ce.gatech.edu/enveomics/index) to build a [recruitment plot](http://enve-omics.ce.gatech.edu/enveomics/docs?t=BlastTab.recplot2.R). There's even a [GUI](http://enve-omics.ce.gatech.edu/enveomics/gui) to make things a bit easier.


## Step 04: Generate summary plots.

```bash
# To Display the program description and parameter options
python 04a_MagicBlast_CoverageMagic_SummaryPlot.py -h

# Example execution:
python 04a_MagicBlast_CoverageMagic_SummaryPlot.py -pre {outfile_prefix} -thd 95 -tad 80
```

Example plot:
![alt text](04b_Example_plot.png "Example summary plot.")


## Step 05: Build Data Table of Whole Genome stats for multiple genomes across multiple metagenomes

Move all the genome.tsv files you want to place in the table into their own directory.

```bash
# To Display the program description and parameter options
python 05_MagicBlast_CoverageMagic_CombineGenomeStats.py -h

# Example execution:
python 05_MagicBlast_CoverageMagic_CombineGenomeStats.py -gtd {genome.tsv directory} -o {outfile_name}
```