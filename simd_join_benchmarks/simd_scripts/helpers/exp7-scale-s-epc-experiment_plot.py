#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py

import commons
import re
import statistics
import subprocess
import csv
import matplotlib.pyplot as plt

filename = "../data/scale-s-output.csv"
mb_of_data = 131072

def plot():
    r_sizes_names = [
        '$R_{size}$ < L2',
        'L2 < $R_{size}$ < L3',
        'L3 < $R_{size}$ < EPC',
        'EPC < $R_{size}$'
    ]
    csvf = open(filename, mode='r')
    csvr = csv.DictReader(csvf)
    all_data = list(csvr)
    r_sizes = sorted(set(map(lambda x:float(x['sizeR']), all_data)))
    s_sizes = sorted(set(map(lambda x:float(x['sizeS']), all_data)))
    algos = sorted(set(map(lambda x:x['alg'], all_data)))
    """
    splitted = [[[y['alg'], y['sizeR'], y['sizeS'], y['throughput']] for y in all_data if y['sizeR'] == str(x)] for x in r_sizes]
    titles = ['R < L2', 'L2 < R < L3', 'L3 < R < EPC', 'EPC < R']
    fig,a = plt.subplots(2,2,figsize=(10,5))
    for i in range(0, len(r_sizes)):
        ds = splitted[i]
        ds = [[[y[0], y[1], y[2], y[3]] for y in ds if y[0] == x] for x in algos]
        x = 1 if i & (1<<1) else 0
        y = 1 if i & (1<<0) else 0
        for j in range(0, len(algos)):
            a[x][y].plot(s_sizes, list(map(lambda x: float(x[3]),ds[j])),
                         '-o', label=algos[j], color=commons.color_alg(algos[j]))
        a[x][y].legend()
        a[x][y].set_xlabel('S size [MB]')
        a[x][y].set_ylabel('Throughput [M rec/s]')
        a[x][y].set_title(titles[i] + ' (' + str(r_sizes[i]) + ' MB)')
        a[x][y].set_ylim([0,110])
    commons.savefig('img/scale-s.png')

    # print graphs per algorithm
    fig = plt.figure(figsize=(8,6))
    plt.clf()
    for alg in algos:
        data = list(filter(lambda x: x['alg'] == alg, all_data))
        data_splitted = [[y for y in data if y['sizeR'] == str(x)] for x in r_sizes]
        plt.subplot(2,2,algos.index(alg)+1)
        for i in range(0, len(r_sizes)):
            plt.plot(s_sizes, list(map(lambda x: float(x['throughput']), data_splitted[i])),
                     '-o', label=r_sizes_names[i], color=commons.color_size(i))
        if alg == 'RSM':
            plt.legend(fontsize='small')
        plt.xlabel('S size [MB]')
        plt.ylabel('Throughput [M rec/s]')
        plt.title(alg)
        plt.ylim([0,110])
        plt.gca().yaxis.grid(linestyle='dashed')
    commons.savefig('img/scale-s-algos.png')
"""
    
    # print only CHT
    # plt.clf()
    fig = plt.figure(figsize=(4,3))
    data = list(filter(lambda x: x['alg'] == 'CHT', all_data))
    data_splitted = [[y for y in data if y['sizeR'] == str(x)] for x in r_sizes]
    markers = ['o', 'v', 'D', 's']
    for i in range(0, len(r_sizes)):
        x_sizes = list(filter(lambda x: x['alg'] == 'CHT', all_data))
        x_sizes = sorted(set(map(lambda x:float(x['sizeS']), x_sizes)))
        y = list(map(lambda x: float(x['throughput']), data_splitted[i]))
        plt.plot(x_sizes, y,
                 label=r_sizes_names[i], color=commons.color_size(i), linewidth=2,
                 marker=markers[i], markersize=8, markeredgecolor='black',
                 markeredgewidth=0.5)
    # plt.legend(fontsize='x-small')
    lines, labels = fig.axes[-1].get_legend_handles_labels()
    fig.legend(lines, labels, fontsize=10, frameon=0,
               ncol=2, bbox_to_anchor = (0.05, 0.95), loc='lower left', borderaxespad=0)
    plt.gca().yaxis.grid(linestyle='dashed')

    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('Size of inner table [MB]',fontsize=10)
    plt.ylabel('Throughput [M rec/s]',fontsize=10)
    plt.title('CHT',fontsize = 12)
    # plt.xlim(left=0)
    plt.ylim([0,2000])
    plt.tight_layout()
    commons.savefig('../img/Fig-e7-scale-s-CHT.pdf')
    
    
    fig = plt.figure(figsize=(4,3))
    data = list(filter(lambda x: x['alg'] == 'CHTfsimd', all_data))
    data_splitted = [[y for y in data if y['sizeR'] == str(x)] for x in r_sizes]
    markers = ['o', 'v', 'D', 's']
    for i in range(0, len(r_sizes)):
        x_sizes = list(filter(lambda x: x['alg'] == 'CHTfsimd', all_data))
        x_sizes = sorted(set(map(lambda x:float(x['sizeS']), x_sizes)))
        y = list(map(lambda x: float(x['throughput']), data_splitted[i]))
        plt.plot(x_sizes, y,
                 label=r_sizes_names[i], color=commons.color_size(i), linewidth=2,
                 marker=markers[i], markersize=8, markeredgecolor='black',
                 markeredgewidth=0.5)
    # plt.legend(fontsize='x-small')
    lines, labels = fig.axes[-1].get_legend_handles_labels()
    fig.legend(lines, labels, fontsize=10, frameon=0,
               ncol=2, bbox_to_anchor = (0.05, 0.95), loc='lower left', borderaxespad=0)
    plt.gca().yaxis.grid(linestyle='dashed')

    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('Size of inner table [MB]',fontsize=10)
    plt.ylabel('Throughput [M rec/s]',fontsize=10)
    plt.title('CHTfsimd', fontsize = 12)
    # plt.xlim(left=0)
    plt.ylim([0,2000])
    plt.tight_layout()
    commons.savefig('../img/Fig-e7-scale-s-CHTfsimd.pdf')

if __name__ == '__main__':
    plot()