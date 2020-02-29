#!/usr/bin/env python

''' Renames sequences in a fastq formatted file.

Magic Blast cuts the read name at the first space character and
reports the name before the space as the query ID.

However, fastq files can be named as:
    @D00468:261:HYTMHBCX2:1:1101:9119:31637 1:N:0:CAGAGAGG+ACTGCATA
    @D00468:261:HYTMHBCX2:1:1101:9119:31637 2:N:0:CAGAGAGG+ACTGCATA
where the unique identifier can come after the space chacter.

Filtering for best hits and retrieving the nucleotide sequence of a magic
blast match becomes much more involved once the unique identifier is lost.

This script renames the fastq sequences with a simple name:
@Read_0000000000001
@Read_0000000000002
etc...

!!! CAUTION !!! Renames file in place (ie: overwrites the input file) !!!

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 7th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys, subprocess

def read_fastq(fq):
    #linecount = 0
    while fq:
        # fastq format line 1x is read name
        name = fq.readline().rstrip()
        # Check for end of file and break.
        if not name: break
        # fastq format line 2x is nucleotide sequence
        seq = fq.readline().rstrip()
        # fastq format line 3x is wierd blank '+'' line
        blank = fq.readline().rstrip()
        # fastq format line 4x is the quality score line
        qal = fq.readline().rstrip()
        # Return these four lines each iteration
        yield (name, seq, blank, qal)

        #linecount += 4
        #if linecount % 4 == 0: yield (name, seq, blank, qal)

def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--fastq_input_file',
        help='Please specify the fastq file to read!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--prefix',
        help='Please specify the prefix to rename with!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    infile = args['fastq_input_file']
    # define temporary file name.
    outfile = args['fastq_input_file'].split('.')[0] + '.temp'
    prefix = args['prefix']

    # open the file to check fastq format
    with open(infile, 'r') as file:
        # only need to read first line
        line = file.readline()
        # Verify fastq format or exit with message
        if not line.startswith('@'):
            print('\n\nError reading fastq formatted read names!\n\n')
            sys.exit('\n\nAre you sure this is a fastq formatted file?\n\n')

    # initialize read count for new read names.
    readcount = 0
    # Read through the input fastq file, change read names, and write new file.
    with open(infile, 'r') as f, open (outfile, 'w') as o:
        # Iterate through fastq format with read_fastq function
        for name, seq, blank, qal in read_fastq(f):
            # increment readcount and write file with numbered read names.
            readcount += 1
            o.write(f'@{prefix}_{readcount:013}\n{seq}\n{blank}\n{qal}\n')

    # Move the temporary file to overwrite the orignal file.
    _ = subprocess.run(['mv', outfile, infile])

if __name__ == "__main__":
    main()
