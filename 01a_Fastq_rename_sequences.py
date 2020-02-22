#!/usr/bin/env python

''' Renames sequences in a fastq file

Magic Blast cuts the read name at the first space character and
reports the name before the space as the query ID.

However, Fastq files can be named as:
    @D00468:261:HYTMHBCX2:1:1101:9119:31637 1:N:0:CAGAGAGG+ACTGCATA
    @D00468:261:HYTMHBCX2:1:1101:9119:31637 2:N:0:CAGAGAGG+ACTGCATA
where the unique identifier can come after the space chacter.

Filtering for best hits and retrieving the fasta sequence of a magic
blast match becomes impossible once the unique identifier is lost.

This script renames the fastq sequences with a simple name:
@Read_0000000000001
@Read_0000000000002
etc...

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

import argparse, subprocess

def read_fastq(fq):
    linecount = 0
    while fq:
        name = fq.readline().rstrip().split(' ')[0][1:]
        if not name: break
        seq = fq.readline().rstrip()
        blank = fq.readline().rstrip()
        qal = fq.readline().rstrip()
        linecount += 4
        if linecount % 4 == 0: yield (name, seq, blank, qal)

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
    outfile = args['fastq_input_file'].split('.')[0] + '.temp'
    prefix = args['prefix']
    readcount = 0

    with open(infile, 'r') as f, open (outfile, 'w') as o:

        for name, seq, blank, qal in read_fastq(f):
            readcount += 1
            o.write(f'@{prefix}_{readcount:013}\n{seq}\n{blank}\n{qal}\n')

    _ = subprocess.run(['mv', outfile, infile])

if __name__ == "__main__":
    main()
