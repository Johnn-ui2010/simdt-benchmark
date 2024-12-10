#Modified the original script exp8-phases-experiment.py
#!/usr/bin/python3
import getopt
import subprocess
import re
import sys

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import commons
import statistics
import yaml


fname_phases = "../data/phases-runtime-output.csv"

"""From the top to down: For CHT (Probe, build, partition)
    For RHT,PRHO_avx (Probe+build, partition)
"""
def plot_phases_per_cycle():
    plot_fname = "../img/Fig-e9-CPU-cycles-per-tuple"
    csvf = open(fname_phases, mode='r')
    csvr = csv.DictReader(csvf)
    all_data = list(csvr)
    all_data = list(filter(lambda x: x['alg'] != 'PRHO', all_data)) #Filter out PRHO
    algos = sorted(set(map(lambda x:x['alg'], all_data)))
    modes = sorted(set(map(lambda x:x['mode'], all_data)))
    all_data = list(filter(lambda x:x['phase'] != "Total", all_data))
    # Delete manually "Phase Join" for CHT and CHT join
    datasets = commons.get_test_dataset_names()

    # a graph for both native and sgx
    # to_modes = [[y for y in all_data if y['mode'] == x] for x in modes]
    # plt.figure(figsize=(10,10))
    # plt.clf()
    # for m in range(0, len(modes)):
    #     mode = modes[m]
    #     to_datasets = [[y for y in to_modes[m] if y['ds'] == x] for x in datasets]
    #     for d in range(0, len(datasets)):
    #         ds = datasets[d]
    #         to_algos = [[y for y in to_datasets[d] if y['alg'] == x] for x in algos]
    #         plt.subplot(2,2,2*m+d+1)
    #         i = 0
    #         agg = 0
    #         for a in range(0, len(algos)):
    #             for j in range(0, len(to_algos[a])):
    #                 numtuples = 6553600 if ds == 'cache-fit' else 65536000
    #                 val = float(to_algos[a][j]['cycles']) / numtuples
    #                 alpha = 1 - 0.7*(j%2)
    #                 plt.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
    #                         color=commons.color_alg(algos[a]), alpha=alpha)
    #                 agg = agg + val
    #             agg = 0
    #             i = i+1
    #
    #         # plt.xlabel("Join algorithm")
    #         plt.ylabel("CPU cycles / tuple")
    #         plt.xticks(np.arange(len(algos)),algos)
    #         plt.title(mode + " - dataset " + ds)
    #         plt.tight_layout()
    # commons.savefig(plot_fname + ".png")

    # now just SGX
    plt.rc('axes', axisbelow=True)
    plt.rcParams.update({'font.size': 15})
    fig = plt.figure(figsize=(8,3.5))
    outer = gridspec.GridSpec(1,2)
    # first dataset
    d = 0
    ds = datasets[d]
    ax = plt.Subplot(fig, outer[d])
    to_datasets = [[y for y in all_data if y['ds'] == x] for x in datasets]
    to_algos = [[y for y in to_datasets[d] if y['alg'] == x] for x in algos]
    i = 0
    agg = 0
    ax.yaxis.grid(linestyle='dashed')
    for a in range(len(algos)):
        alg = algos[a]
        for j in range(len(to_algos[a])):
            if ds == 'cache-fit':
                numtuples = 6553600
            elif ds == 'cache-exceed':
                numtuples = 65536000
            else:
                raise ValueError('Unknown dataset: ' + ds)

            val = float(to_algos[a][j]['cycles']) / numtuples
            alpha = 1 - 0.7*(j%2)
            # first print a white background to cover the grid lines
            ax.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
                   color='white', edgecolor='black')
            #now print the real bar
            ax.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
                   color= commons.color_alg(alg), alpha=alpha, edgecolor='black')
            agg = agg + val
        agg = 0
        i = i+1
    ax.set_xticks(np.arange(len(algos)))
    ax.set_xticklabels(algos, rotation=45, fontsize=9)
    # ax.set_xlabel("Join algorithm")
    ax.set_ylabel("CPU cycles / tuple")
    ax.set_ylim([0, 150])
    ax.set_title('(' + chr(97+d) + ") Dataset $\it{" + ds + "}$", y=-0.4)
    #ax.legend()
    fig.add_subplot(ax)

    #second dataset
    d = 1
    ds = datasets[d]
    to_algos = [[y for y in to_datasets[d] if y['alg'] == x] for x in algos]
    ax2 = plt.Subplot(fig, outer[1])
    i = 0
    agg = 0
    for a in range(len(algos)):
        alg = algos[a]
        for j in range(0, len(to_algos[a])):
            if ds == 'cache-fit':
                numtuples = 6553600
            elif ds == 'cache-exceed':
                numtuples = 65536000
            else:
                raise ValueError('Unknown dataset: ' + ds)
            val = float(to_algos[a][j]['cycles']) / numtuples
            alpha = 1 - 0.7*(j%2)
            # first print a white background to cover the grid lines
            #ax.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
             #      color='white', edgecolor='black')
            ax2.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
                    color='white', edgecolor='black')
            #now print the real bar
            #ax.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
             #      color=commons.color_alg(alg), alpha=alpha, edgecolor='black')
            ax2.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
                    color=commons.color_alg(alg), alpha=alpha, edgecolor='black')
            agg = agg + val
            
        agg = 0
        i = i+1

    ax2.set_ylim(0, 150)
    #ax.set_ylim(500, 15500)
    #ax.set_yticks((1000,5000,10000,15000))
    # hide the spines between ax and ax2
    #ax.spines['bottom'].set_visible(False)
    #ax2.spines['top'].set_visible(False)

    #ax.xaxis.tick_top()
    #ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()

    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    #kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)

    #kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    #ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    #ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    ax2.yaxis.grid(linestyle='dashed')
    ax2.set_xticks(np.arange(len(algos)))
    ax2.set_xticklabels(algos, rotation=45, fontsize=9)

    #ax2.legend()
    #ax2.set_xlabel("Join algorithm")
    ax2.set_title('(' + chr(97+1) + ") Dataset $\it{" + ds + "}$", y=-0.4)
    #fig.text(0.5, 0.57, 'CPU cycles [B]', va='center', rotation='vertical')

    #fig.add_subplot(ax)
    fig.add_subplot(ax2)

    plt.tight_layout()
    commons.savefig(plot_fname + '.pdf')


if __name__ == '__main__':

    plot_phases_per_cycle()
