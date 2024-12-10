#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# filename_more_r: R:S tuples = 1:4
# filename_more_s: R:S tuples = 4:1

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

#Not finished.
filename_more_r_phases = "../data/CHT_more_r_phases.csv"
#filename_more_r_phases = "../data/CHT_more_r_phases_small.csv"

filename_more_s_phases = "../data/CHT_more_s_phases.csv"
mb_of_data = 131072

"""From the top to down: For CHT (Probe, build, partition)
"""
def plot_phases_per_cycle():
    plot_fname = "../img/Fig-e11-CHT_more_r_phases"
    csvf = open(filename_more_r_phases, mode='r')
    csvr = csv.DictReader(csvf)
    all_data = list(csvr)
    #algos = sorted(set(map(lambda x:x['alg'], all_data)))
    algos = ['CHT', 'CHTfsimd']
    print(algos)
    modes = sorted(set(map(lambda x:x['mode'], all_data)))
    all_data = list(filter(lambda x:x['phase'] != "Total", all_data))
    # Delete manually "Phase Join" for CHT and CHT join

    # a graph for both native and sgx
    """
    to_modes = [[y for y in all_data if y['mode'] == x] for x in modes]
    plt.figure(figsize=(5,5))
    plt.clf()
    for m in range(0, len(modes)):
        mode = modes[m]
        to_datasets = [y for y in to_modes[m]]
        to_algos = [[y for y in to_datasets if y['alg'] == x] for x in algos]
        plt.subplot(2,2,2*m+1)
        i = 0
        agg = 0
        for a in range(0, len(algos)):
            for j in range(0, len(to_algos[a])):
                numtuples = 52829996
                val = float(to_algos[a][j]['cycles']) / numtuples
                alpha = 1 - 0.7*(j%2)
                plt.bar(i, val, bottom=agg, label=to_algos[a][j]['phase'],
                        color=commons.color_alg(algos[a]), alpha=alpha)
            agg = agg + val
            agg = 0
            i = i+1
            # plt.xlabel("Join algorithm")
        plt.ylabel("CPU cycles / tuple")
        plt.xticks(np.arange(len(algos)),algos)
        plt.title(mode + " - real dataset ")
        plt.tight_layout()
    commons.savefig(plot_fname + ".pdf")
    """
    # now just SGX
    plt.rc('axes', axisbelow=True)
    plt.rcParams.update({'font.size': 15})
    fig = plt.figure(figsize=(8,3.5))
    outer = gridspec.GridSpec(1,2)
    # first dataset
    d = 0
    ax = plt.Subplot(fig, outer[d])
    to_datasets = all_data
    to_datasets = [y for y in all_data if y['mode'] == "sgx"]
    to_algos = [[y for y in to_datasets if y['alg'] == x] for x in algos]
    i = 0
    agg = 0
    ax.yaxis.grid(linestyle='dashed')
    numtuples = 52829996

    for a in range(len(algos)):
        alg = algos[a]
        for j in range(len(to_algos[a])):
            print(to_algos[a][j])
            print(to_algos[a][j]['cycles'])
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
    ax.set_ylim([0, 75])
    #ax.set_title('(' + chr(97+d) + ") dataset", y=-0.4)
    #ax.legend()
    fig.add_subplot(ax)

    plt.tight_layout()
    commons.savefig(plot_fname + '.pdf')

if __name__ == '__main__':

    plot_phases_per_cycle()