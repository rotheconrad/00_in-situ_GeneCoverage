#!/usr/bin/env python

'''Calculate and plot histogram of pIdent from tabular Magic Blast.

This tool takes the following input parameters:

    * uniqueID.blast - tabular Magic Blast File

This script returns the following files:

    * Histogram as tsv file
    * Histogram Plot .png format

This script requires the following packages:

    * argparse
    * matplotlib

This file can also be imported as a module and contains the follwing 
functions:

    * parse_magic_blast - parses tabular magic blast file
    * calc_pIdent_hist - calculates the histogram
    * parse_hist_file - parses the histogram data
    * plot_pIdent_hist - plots the data
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Wednesday, August 28th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def parse_magic_blast(file, data_dict):
    """ parses file returns updated dictionary data_dict """

    alignment_lengths = {i: 0 for i in range(70,101)}
    query_lengths = {i: 0 for i in range(70,101)}
    read_counts = {i: 0 for i in range(70,101)}

    name = file.split('_')[0]

    with open(file, 'r') as f:
        for l in f:
            if l.startswith('#'): continue
            X = l.rstrip().split('\t')
            pident = int(X[2].split('.')[0])
            astrt = min(int(X[8]), int(X[9]))
            astp = max(int(X[8]), int(X[9]))
            aLen = astp - astrt # read alignment length
            qLen = int(X[15]) # full length of read

            if pident >= 70:
                alignment_lengths[pident] += aLen
                query_lengths[pident] += qLen
                read_counts[pident] += 1

        data_dict['alen'] = alignment_lengths
        data_dict['qlen'] = query_lengths
        data_dict['rcount'] = read_counts

    return data_dict


def calc_pIdent_hist(infile, outfile):
    """Reads the files, orchestrates the parsing, calculation and output.

    Parameters
    ----------
    input file : str
        Tabular Blast output file names separated by spaces
    out file : str
        Prefix to use when naming the output files

    Returns
    -------
    Tab delimimted data table file.
        
    """
    data_dict = {}

    print(f'Parsing file {infile}')
    parse_magic_blast(infile, data_dict)


    with open(outfile, 'w') as o:

        print(f'File parsed. Writing to {outfile}')
        header_order = sorted(data_dict.keys())
        header = 'pIdent\t' + '\t'.join(header_order) + '\n'

        o.write(header)

        for i in reversed(range(70, 101)):
            buildLine = [str(i)]

            for j in header_order:
                buildLine.append(str(data_dict[j][i]))

            o.write('\t'.join(buildLine) + '\n')


def parse_hist_file(infile):
    """ Parses the file into arrays to be plotted """

    d = {'ys': [], 'xs': []}

    with open(infile, 'r') as f:
        _ = f.readline() #skip header
        for l in f:
            X = l.split('\t')
            y = X[0]
            x = int(X[1])
            d['ys'].append(y)
            d['xs'].append(x)

    return d


def plot_pIdent_hist(d, outfile):
    """ Takes d from parse_files function and plots the data """

    print('Plotting...')
    # Set the colors
    bar_color = '#2171b5'
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    alpha = 0.6


    # Build the plot
    fig, ax = plt.subplots(figsize=(7, 10))

    '''
    # Plot titles
    ax.set_title(
        f'Histogram of Percent Identity\nof Read Alignments',
        fontsize=20, y=1.02
        )
    '''

    # Plot labels
    ax.set_xlabel(
                'Number of Base Pairs Aligned',
                fontsize=14, fontweight='bold'
                )
    ax.set_ylabel(
                'Percent Identity of Alignment',
                fontsize=14, fontweight='bold'
                )

    # Set plot/grid style
    ax.minorticks_on()
    ax.tick_params(
        which='minor', axis='both', left=False, bottom=False
        )
    ax.tick_params(
                which='major', axis='both',
                left=False, bottom=True,
                size=4, width=2, tickdir='in',
                labelsize=11, zorder=10
                )
    ax.xaxis.grid(
        which="minor", color=gridm, linestyle='--',
        linewidth=1, alpha=0.6, zorder=1
        )
    ax.xaxis.grid(
        which="major", color=gridM, linestyle='--',
        linewidth=1.5, alpha=0.4, zorder=1
        )
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_linewidth(2)

    # Plot the data
    ax.barh(
        d['ys'],
        d['xs'],
        align='center',
        height=0.9,
        color=bar_color,
        alpha=alpha,
        )

    # Set plot axis ranges
    #ax.set_xlim(left=0, right=int((max(d['xs'])+min(d['xs']))))

    # adjust layout, save, and close
    plt.gca().invert_yaxis()
    fig.set_tight_layout(True)
    plt.savefig(outfile)
    plt.close()

    print('\n\nComplete success space cowboy! Hold on to your boots.\n\n')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the pIdent_hist.tsv file!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\nRunning Script...\n')
    outpre = args['input_file'].split('.')[0]
    out_hist = f"{outpre}.hist"
    _ = calc_pIdent_hist(args['input_file'], out_hist)
    d = parse_hist_file(out_hist)
    out_plot = f"{outpre}.hist.png"
    _ = plot_pIdent_hist(d, out_plot)


if __name__ == "__main__":
    main()
