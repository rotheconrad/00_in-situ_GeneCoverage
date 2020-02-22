#!/usr/bin/env python

''' Rename fasta sequences

Reads fasta file and replaces sequence names with >prefix_i where 
i = 1 through the number of sequences in the fasta file.

Renames file in place (ie: overwrites the input file).

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

import argparse, subprocess

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

def Fasta_rename_sequences(infile, prefix):

    outfile = prefix + '.temp'

    with open(infile, 'r+') as f, open(outfile, 'w') as o:
        i = 1
        for name, seq in read_fasta(f):
            newName = f'>{prefix}_{i:013}\n{seq}\n'
            o.write(newName)
            i += 1

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

    # Run the renaming function:
    # Print renaming file
    Fasta_rename_sequences(args['input_file'], args['naming_prefix'])
    
if __name__ == "__main__":
    main()

