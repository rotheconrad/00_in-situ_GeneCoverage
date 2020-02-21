#!/usr/bin/env python

'''Plot Results from MagicBlast02b_CoveragePlus.py

This tool takes the following input parameters:

    * The -o out_file_prefix from CoveragePlus script used to name:
        *_genome_by_bp.tsv
        *_contig_ani.tsv
        *_contig_tad.tsv
        *_gene_ani.tsv
        *_gene_tad.tsv

This script returns the following files:

    * A .png file of plots for TAD (coverage) and ANIr by genome, contig,
      and genes along with plots for contig and gene length verse TAD.

This script requires the following packages:

    * argparse
    * collections.defaultdict
    * matplotlib
    * scipy.stats.stats.pearsonr

This file can also be imported as a module and contains the follwing 
functions:

    * plot_ANIr_TAD - short description
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, August 15th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr


def parse_files(gbp, cani, ctad, gani, gtad):
    """ Parses the files into arrays to be plotted """

    gene_xs_dict = defaultdict(int)

    d = {
            'gbp_xs': [],
            'gbp_ys': [],
            'tad_cxs': 0,
            'tad_gxs': 0,
            'ani_cxs': 0,
            'ani_gxs': 0,
            'wcanis': [],
            'wganis': [],
            'wctads': [],
            'wgtads': [],
            'canis': [],
            'ganis': [],
            'ctads': [],
            'gtads': [],
            'aniclens': [],
            'tadclens': [],
            'aniglens': [],
            'tadglens': []
        }

    print(f'Reading {gbp}...')
    with open(gbp, 'r') as f:
        _ = f.readline() # skip header
        for l in f:
            X = l.rstrip().split('\t')
            d['gbp_xs'].append(int(X[0]))
            d['gbp_ys'].append(int(X[1]))

    print(f'Reading {cani}...')
    with open(cani, 'r') as f:
        _ = f.readline() # skip header
        for i, l in enumerate(f):
            X = l.rstrip().split('\t')
            value = float(X[1])
            length = int(X[2])
            d['ani_cxs'] += length
            d['wcanis'].extend([value]*length)
            d['canis'].append(value)
            d['aniclens'].append(length)

    print(f'Reading {ctad}...')
    with open(ctad, 'r') as f:
        _ = f.readline() # skip header
        for i, l in enumerate(f):
            X = l.rstrip().split('\t')
            value = float(X[1])
            length = int(X[2])
            d['tad_cxs'] += length
            d['wctads'].extend([value]*length)
            d['ctads'].append(value)
            d['tadclens'].append(length)

    print(f'Reading {gani}...')
    with open(gani, 'r') as f:
        _ = f.readline() # skip header
        for i, l in enumerate(f):
            X = l.rstrip().split('\t')
            contig_name = int(X[0].split('_')[-2])
            value = float(X[1])
            length = int(X[2])
            d['ani_gxs'] += length
            d['wganis'].extend([value]*length)
            d['ganis'].append(value)
            d['aniglens'].append(length)

    print(f'Reading {gtad}...')
    with open(gtad, 'r') as f:
        _ = f.readline() # skip header
        for i, l in enumerate(f):
            X = l.rstrip().split('\t')
            value = float(X[1])
            length = int(X[2])
            d['tad_gxs'] += length
            d['wgtads'].extend([value]*length)
            d['gtads'].append(value)
            d['tadglens'].append(length)

    return d


def plot_ANIr_TAD(d, tad, thd, outfile):
    """ Takes d from parse_files function and plots the data """

    # Caculate Correlations
    print('\nComputing correlations...\n')
    v_cani_corr = pearsonr(d['aniclens'], d['canis'])
    v_gani_corr = pearsonr(d['aniglens'], d['ganis'])
    v_ctad_corr = pearsonr(d['tadclens'], d['ctads'])
    v_gtad_corr = pearsonr(d['tadglens'], d['gtads'])

    # Correlation Strings to Print
    cani_corr = (
        f'Pearson r: {round(v_cani_corr[0], 3)}, '
        f'p={round(v_cani_corr[1], 3)}'
        )
    gani_corr = (
        f'Pearson r: {round(v_gani_corr[0], 3)}, '
        f'p={round(v_gani_corr[1], 3)}'
        )
    ctad_corr = (
        f'Pearson r: {round(v_ctad_corr[0], 3)}, '
        f'p={round(v_ctad_corr[1], 3)}'
        )
    gtad_corr = (
        f'Pearson r: {round(v_gtad_corr[0], 3)}, '
        f'p={round(v_gtad_corr[1], 3)}'
        )

    print('Plotting...')
    # Set the colors
    gbp_color = '#4d4d4d'
    cani_color = '#af8dc3'
    gani_color = '#af8dc3'
    ctad_color = '#7fbf7b'
    gtad_color = '#7fbf7b'
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    alpha = 0.25
    alpha2 = 0.75

    # Build the plot
    fig, (
        (ax1, ax2),
        (ax3, ax4),
        (ax5, ax6),
        (ax7, ax8),
        (ax9, ax10)
        ) = plt.subplots(
                        5, 2,
                        figsize=(25, 15),
                        gridspec_kw={'width_ratios': [5, 1]},
                        )
    ax1.get_shared_y_axes().join(ax1, ax3, ax5)
    ax7.get_shared_y_axes().join(ax7, ax8, ax9, ax10)

    # Plot titles
    ax1.set_title(
        f'Sequencing Depth by Base Pairs (bps)',
        fontsize=20, y=1.02
        )
    ax2.set_title(
        f'Per Base Pair Sequencing Depth',
        fontsize=20, y=1.02
        )
    ax3.set_title(
        f'TAD_{tad} by Contig',
        fontsize=20, y=1.02
        )
    ax4.set_title(
        f'TAD by Contig Length',
        fontsize=20, y=1.02
        )
    ax5.set_title(
        f'TAD_{tad} by Gene',
        fontsize=20, y=1.02
        )
    ax6.set_title(
        f'TAD by Gene Length',
        fontsize=20, y=1.02
        )
    ax7.set_title(
        f'ANIr_{thd} by Contig',
        fontsize=20, y=1.02
        )
    ax8.set_title(
        f'ANIr by Contig Length',
        fontsize=20, y=1.02
        )
    ax9.set_title(
        f'ANIr_{thd} by Gene',
        fontsize=20, y=1.02
        )
    ax10.set_title(
        f'ANIr by Gene Length',
        fontsize=20, y=1.02
        )

    # Plot labels
    ax9.set_xlabel('Genome Position (bps)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Sequencing Depth', fontsize=14, fontweight='bold')
    #ax3.set_xlabel('Contig Number', fontsize=14, fontweight='bold')
    #ax5.set_xlabel('Gene Number', fontsize=14, fontweight='bold')
    #ax7.set_xlabel('Contig Number', fontsize=14, fontweight='bold')
    #ax9.set_xlabel('Gene Number', fontsize=14, fontweight='bold')
    ax10.set_xlabel('Length in Base Pairs', fontsize=14, fontweight='bold')

    ax1.set_ylabel('TAD', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Count', fontsize=14, fontweight='bold')
    ax3.set_ylabel('TAD', fontsize=14, fontweight='bold')
    ax5.set_ylabel('TAD', fontsize=14, fontweight='bold')
    ax7.set_ylabel('ANIr', fontsize=14, fontweight='bold')
    ax9.set_ylabel('ANIr', fontsize=14, fontweight='bold')

    # Correlation Text
    ax4.text(
        0.99, 0.96,
        cani_corr,
        verticalalignment='top',
        horizontalalignment='right',
        transform=ax4.transAxes,
        fontsize=12
        )
    ax6.text(
        0.99, 0.96,
        gani_corr,
        verticalalignment='top',
        horizontalalignment='right',
        transform=ax6.transAxes,
        fontsize=12
        )
    ax8.text(
        0.99, 0.04,
        ctad_corr,
        verticalalignment='bottom',
        horizontalalignment='right',
        transform=ax8.transAxes,
        fontsize=12
        )
    ax10.text(
        0.99, 0.04,
        gtad_corr,
        verticalalignment='bottom',
        horizontalalignment='right',
        transform=ax10.transAxes,
        fontsize=12
        )

    # Set plot/grid style
    for ax in fig.axes:
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
        ax.yaxis.grid(
            which="minor", color=gridm, linestyle='--',
            linewidth=1, alpha=0.6, zorder=1
            )
        ax.yaxis.grid(
            which="major", color=gridM, linestyle='--',
            linewidth=1.5, alpha=0.4, zorder=1
            )
        ax.set_axisbelow(True)
        for spine in ax.spines.values(): spine.set_linewidth(2)

    # Plot the data
    print(f"... Sequencing Depth by {len(d['gbp_xs']):,} base pair...")
    ax1.plot(
        d['gbp_xs'],
        d['gbp_ys'],
        color=gbp_color,
        linestyle='-',
        lw=0.1,
        alpha=alpha2,
        label='Sequencing Depth'
        )
    bins = int(max(d['gbp_ys']))
    print('... Histogram of per base pair sequencing depths...')
    print(
        f"... ... {len(d['gbp_ys']):,} values to organize into {bins} bins..."
        )
    ax2.hist(
        d['gbp_ys'],
        bins,
        color=gbp_color,
        alpha=alpha2
        )
    print('... TAD and ANI by contig & gene...')
    ax3.plot(
        range(0, d['tad_cxs']),
        d['wctads'],
        color=ctad_color,
        linestyle='-',
        lw=1,
        label='Contig TAD'
        )
    ax4.scatter(
        d['tadclens'],
        d['ctads'],
        marker='o',
        s=20,
        color=ctad_color,
        alpha=alpha
        )
    ax5.plot(
        range(0, d['tad_gxs']),
        d['wgtads'],
        color=gtad_color,
        linestyle='-',
        lw=0.5,
        label='Gene TAD'
        )
    ax6.scatter(
        d['tadglens'],
        d['gtads'],
        marker='s',
        s=20,
        color=gtad_color,
        alpha=alpha
        )
    ax7.plot(
        range(0, d['ani_cxs']),
        d['wcanis'],
        color=cani_color,
        linestyle='-',
        lw=1,
        label='Contig ANIr'
        )
    ax8.scatter(
        d['aniclens'],
        d['canis'],
        marker='o',
        s=20,
        color=cani_color,
        alpha=alpha,
        )
    ax9.plot(
        range(0, d['ani_gxs']),
        d['wganis'],
        color=gani_color,
        linestyle='-',
        lw=0.5,
        label='Gene ANIr'
        )
    ax10.scatter(
        d['aniglens'],
        d['ganis'],
        marker='s',
        s=20,
        color=gani_color,
        alpha=0.5,
        )

    # Set plot axis ranges
    ax1.set_xlim(left=-5, right=max(d['gbp_xs'])+5)
    ax2.set_xlim(left=-5, right=max(d['gbp_ys'])+5)
    ax3.set_xlim(left=-5, right=d['tad_cxs']+5)
    ax4.set_xlim(left=-5, right=max(d['tadclens'])+5)
    ax5.set_xlim(left=-5, right=d['tad_gxs']+5)
    ax6.set_xlim(left=-5, right=max(d['tadglens'])+5)
    ax7.set_xlim(left=-5, right=d['ani_cxs']+5)
    ax8.set_xlim(left=-5, right=max(d['aniclens'])+5)
    ax9.set_xlim(left=-5, right=d['ani_gxs']+5)
    ax10.set_xlim(left=-5, right=max(d['aniglens'])+5)

    ax7.set_ylim(bottom=thd-0.5, top=100.5)

    # adjust layout, save, and close
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
        '-pre', '--input_file_prefix',
        help='Please specify the prefix for TAD ANI tsv files!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-thd', '--pIdent_threshold_cutoff',
        help='Please specify pIdent threshold to use! (ie: 95)',
        metavar=':',
        type=float,
        required=True
        )
    parser.add_argument(
        '-tad', '--truncated_avg_depth_value',
        help='Please specify TAD value! (ie: 80 or 90)',
        metavar=':',
        type=float,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\nRunning Script...\n')

    gbp = f"{args['input_file_prefix']}_genome_by_bp.tsv"
    cani = f"{args['input_file_prefix']}_contig_ani.tsv"
    ctad = f"{args['input_file_prefix']}_contig_tad.tsv"
    gani = f"{args['input_file_prefix']}_gene_ani.tsv"
    gtad = f"{args['input_file_prefix']}_gene_tad.tsv"
    out_file = f"{args['input_file_prefix']}_TAD_ANIr_plot.png"

    d = parse_files(gbp, cani, ctad, gani, gtad)

    _ = plot_ANIr_TAD(
                    d,
                    args['truncated_avg_depth_value'],
                    args['pIdent_threshold_cutoff'],
                    out_file
                    )


if __name__ == "__main__":
    main()
