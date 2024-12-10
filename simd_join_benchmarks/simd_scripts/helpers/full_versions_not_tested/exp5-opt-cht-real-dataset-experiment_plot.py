#!/usr/bin/python3
# Modified the original: exp5-real-dataset-experiment.py
import getopt
import subprocess
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import csv
import commons
import statistics

import yaml
from matplotlib.patches import Patch

#Change the path
fname_output = "../data/real-dataset-output_opt.csv"

def plot_throughput():
    #Change the path
    plot_filename = "../img/Figure-Join-algorithms-throughput-with-IMDb_opt"
    csvf = open(fname_output, mode='r')
    csvr = csv.DictReader(csvf)
    all_data = list(csvr)
    algos = list(set(map(lambda x:x['alg'], all_data)))
    modes = sorted(set(map(lambda x:x['mode'], all_data)))
    width = 0.4
    to_modes = [[y for y in all_data if y['mode'] == x] for x in modes]

    # graph per dataset
    plt.rc('axes', axisbelow=True)
    plt.rcParams.update({'font.size': 15})
    fig = plt.figure(figsize=(5,4))
    plt.clf()
    to_modes = [[y for y in all_data if y['mode'] == x] for x in modes]
    for m in range(0, len(modes)):
        plt.gca().yaxis.grid(linestyle='dashed')
        if m == 0:
            br = np.arange(len(algos))
            br = [x - 0.2 for x in br]
        else:
            br = [x + width for x in br]

        label = modes[m]
        hatch = '\\' if modes[m] == 'native' else ''
        to_modes[m] = sorted(to_modes[m], key=lambda x:x['alg'])
        #colors = list(map(lambda x: commons.color_alg(x['alg']), to_modes[m]))
        throughputs = list(map(lambda x: float(x['throughput']), to_modes[m]))
        plt.bar(br, throughputs,width=width, label=label, hatch=hatch,
                edgecolor='black') #color=colors, 
        for x, y in zip(br, list(map(lambda x: float(x['throughput']), to_modes[m]))):
            if y < 4:
                plt.text(x-0.15, y+6, str(y), rotation=90)
        plt.xlabel("Join algorithm")
        plt.ylabel("Throughput [M rec/s]")
        
        #Adding a constant red line with the value of CHT.
        plt.axhline(y=throughputs[0], color='r', linestyle='-')
        # plt.ylim([0, 220])
        plt.xticks(np.arange(len(to_modes[m])), list(map(lambda x:x['alg'],to_modes[m])),
                   rotation=45, size=7)
        # plt.title("IMDb dataset")
    # plt.gca().yaxis.grid(linestyle='dashed')
    
    ax1 = plt.gca()
    ax1.yaxis.grid(linestyle='dashed')
    ax2 = ax1.twinx()
    # ax2.set_yticks([0])
    ax2.tick_params(axis='y', colors='white')
    ax2.set_ylabel(' ')
    legend_elements = [#Patch(label='Hatches:', alpha=0),
                       Patch(facecolor='white', edgecolor='black',
                             hatch='\\\\', label='plain CPU'),
                       Patch(facecolor='white', edgecolor='black',
                             label='TEE')]
    fig.legend(handles=legend_elements, ncol=3, frameon=False,
               bbox_to_anchor=(0.15,0.91,1,0), loc="lower left",
               handletextpad=0.5)
    # plt.annotate('Hatches:',xy=(170, 850), xycoords='figure pixels')
    commons.savefig(plot_filename + ".png", tight_layout=True)


if __name__ == '__main__':

    plot_throughput()
