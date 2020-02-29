#!/usr/bin/env python

''' Renames sequences in a fasta formatted file.

Magic Blast cuts the read name at the first space character and
reports the name before the space as the query ID.

However, fasta files can be named as:
    >D00468:261:HYTMHBCX2:1:1101:9119:31637 1:N:0:CAGAGAGG+ACTGCATA
    >D00468:261:HYTMHBCX2:1:1101:9119:31637 2:N:0:CAGAGAGG+ACTGCATA
where the unique identifier can come after the space chacter.

Filtering for best hits and retrieving the nucleotide sequence of a magic
blast match becomes much more involved once the unique identifier is lost.

This script renames the fasta sequences with a simple name:
>Read_0000000000001
>Read_0000000000002
etc...

!!! CAUTION !!! Renames file in place (ie: overwrites the input file) !!!

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 21st, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys, subprocess

def read_fasta(fp):
    # initialize name, seq objects
    name, seq = None, []
    # iterate through lines in file
    for line in fp:
        # strip new line characters
        line = line.rstrip()
        # check if line is read name
        if line.startswith(">"):
            # return prior name, seq object at each new name
            if name: yield (name, ''.join(seq))
            # define new name, empty seq object for next sequence
            name, seq = line, []
        # else it is sequence so add to current seq object
        else:
            seq.append(line)
    # this returns the last name, seq pair of the file.
    if name: yield (name, ''.join(seq))

def Fasta_rename_sequences(infile, prefix):

    # define temporary file name.
    outfile = prefix + '.temp'

    # open the infile and temp file: read and write.
    with open(infile, 'r+') as f, open(outfile, 'w') as o:
        # initialize read count for new read names.
        i = 1
        # Iterate through the fasta format with read_fasta function.
        for name, seq in read_fasta(f):
            # Write temp file with numbered read names.
            o.write(f'>{prefix}_{i:013}\n{seq}\n')
            # increment read count
            i += 1

    # Move the temp file to replace the original file.
    _ = subprocess.run(['mv', outfile, infile])

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the fasta file to rename deflines!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--naming_prefix',
        help='Please specify the prefix to use for renaming sequences!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # open the file to check fasta format
    with open(args['input_file'], 'r') as file:
        # Only need to read first line
        line = file.readline()
        # Verify fasta format or exit with message.
        if not line.startswith('>'):
            print('\n\nError reading fasta formatted read names!\n\n')
            sys.exit('\n\nAre you sure this is a fasta formatted file?\n\n')

    # Run the renaming function:
    Fasta_rename_sequences(args['input_file'], args['naming_prefix'])
    
if __name__ == "__main__":
    main()

