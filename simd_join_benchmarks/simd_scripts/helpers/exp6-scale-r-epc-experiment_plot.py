#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# Original script modified: exp6-scale-s-experiment-SUDO.py
import commons
import re
import statistics
import subprocess
import csv
import matplotlib.pyplot as plt

filename = "../data/scale-r-output.csv"
mb_of_data = 131072

def plot():
    s_sizes_names = [
        '$S_{size}$ < L2',
        'L2 < $S_{size}$ < L3',
        'L3 < $S_{size}$ < EPC',
        'EPC < $S_{size}$'
    ]
    csvf = open(filename, mode='r')
    csvr = csv.DictReader(csvf)
    all_data = list(csvr)
    r_sizes = sorted(set(map(lambda x:float(x['sizeR']), all_data)))
    s_sizes = sorted(set(map(lambda x:float(x['sizeS']), all_data)))
    algos = sorted(set(map(lambda x:x['alg'], all_data)))
    """
    splitted = [[[y['alg'], y['sizeR'], y['sizeS'], y['throughput']] for y in all_data if y['sizeR'] == str(x)] for x in r_sizes]
    titles = ['S < L2', 'L2 < S < L3', 'L3 < S < EPC', 'EPC < S']
    fig,a = plt.subplots(2,2,figsize=(10,5))
    
    for i in range(0, len(s_sizes)):
        ds = splitted[i]
        ds = [[[y[0], y[1], y[2], y[3]] for y in ds if y[0] == x] for x in algos]
        x = 1 if i & (1<<1) else 0
        y = 1 if i & (1<<0) else 0
        for j in range(0, len(algos)):
            a[x][y].plot(r_sizes, list(map(lambda x: float(x[3]),ds[j])),
                         '-o', label=algos[j], color=commons.color_alg(algos[j]))
        a[x][y].legend()
        a[x][y].set_xlabel('R size [MB]')
        a[x][y].set_ylabel('Throughput [M rec/s]')
        a[x][y].set_title(titles[i] + ' (' + str(s_sizes[i]) + ' MB)')
        a[x][y].set_ylim([0,110])
    commons.savefig('img/scale-r.png')

    # print graphs per algorithm
    fig = plt.figure(figsize=(8,6))
    plt.clf()
    for alg in algos:
        data = list(filter(lambda x: x['alg'] == alg, all_data))
        data_splitted = [[y for y in data if y['sizeS'] == str(x)] for x in s_sizes]
        plt.subplot(2,2,algos.index(alg)+1)
        for i in range(0, len(s_sizes)):
            plt.plot(r_sizes, list(map(lambda x: float(x['throughput']), data_splitted[i])),
                     '-o', label=s_sizes_names[i], color=commons.color_size(i))
        if alg == 'RSM':
            plt.legend(fontsize='small')
        plt.xlabel('R size [MB]')
        plt.ylabel('Throughput [M rec/s]')
        plt.title(alg)
        plt.ylim([0,110])
        plt.gca().yaxis.grid(linestyle='dashed')
    commons.savefig('img/scale-r-algos.png')
    """
    # print only CHT
    # plt.clf()
    fig = plt.figure(figsize=(4,3))
    data = list(filter(lambda x: x['alg'] == 'CHT', all_data))
    data_splitted = [[y for y in data if y['sizeS'] == str(x)] for x in s_sizes]
    markers = ['o', 'v', 'D', 's']
    for i in range(0, len(s_sizes)):
        x_sizes = list(filter(lambda x: x['alg'] == 'CHT', all_data))
        x_sizes = sorted(set(map(lambda x:float(x['sizeR']), x_sizes)))
        y = list(map(lambda x: float(x['throughput']), data_splitted[i]))
        plt.plot(x_sizes, y,
                 label=s_sizes_names[i], color=commons.color_size(i), linewidth=2,
                 marker=markers[i], markersize=8, markeredgecolor='black',
                 markeredgewidth=0.5)
    # plt.legend(fontsize='x-small')
    lines, labels = fig.axes[-1].get_legend_handles_labels()
    fig.legend(lines, labels, fontsize=10, frameon=0,
               ncol=2, bbox_to_anchor = (0.05, 0.95), loc='lower left', borderaxespad=0)
    plt.gca().yaxis.grid(linestyle='dashed')
    
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('Size of outer table [MB]',fontsize=10)
    plt.ylabel('Throughput [M rec/s]',fontsize=10)
    plt.title('CHT', fontsize = 12)
    # plt.xlim(left=0)
    plt.ylim([0,300])
    commons.savefig('../img/Fig-e6-scale-r-CHT.pdf')
    
    fig.clear()
    fig = plt.figure(figsize=(4,3))
    data = list(filter(lambda x: x['alg'] == 'CHTfsimd', all_data))
    data_splitted = [[y for y in data if y['sizeS'] == str(x)] for x in s_sizes]
    markers = ['o', 'v', 'D', 's']
    for i in range(0, len(s_sizes)):
        x_sizes = list(filter(lambda x: x['alg'] == 'CHTfsimd', all_data))
        x_sizes = sorted(set(map(lambda x:float(x['sizeR']), x_sizes)))
        y = list(map(lambda x: float(x['throughput']), data_splitted[i]))
        plt.plot(x_sizes, y,
                 label=s_sizes_names[i], color=commons.color_size(i), linewidth=2,
                 marker=markers[i], markersize=8, markeredgecolor='black',
                 markeredgewidth=0.5)
    # plt.legend(fontsize='x-small')
    lines, labels = fig.axes[-1].get_legend_handles_labels()
    fig.legend(lines, labels, fontsize=10, frameon=0,
               ncol=2, bbox_to_anchor = (0.05, 0.95), loc='lower left', borderaxespad=0)
    plt.gca().yaxis.grid(linestyle='dashed')

    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('Size of outer table [MB]',fontsize=10)
    plt.ylabel('Throughput [M rec/s]',fontsize=10)
    plt.title('CHTfsimd',fontsize = 12)
    # plt.xlim(left=0)
    plt.ylim([0,300])
    commons.savefig('../img/Fig-e6-scale-r-CHTfsimd.pdf')

if __name__ == '__main__':
    plot()
