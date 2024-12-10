#!/usr/bin/python3
#  Original script modified: exp1-off-the-shelf-performance.py

import getopt
import subprocess
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import csv
import commons
import statistics
from matplotlib.patches import Patch
from collections import OrderedDict
import yaml

fname_throughput = "../data/throughput-output.csv"


def run_join(mode, alg, ds, threads, reps):

    f_throughput = open(fname_throughput, 'a')
    throughput_array = []
    throughput = ''
    dic_phases = {}

    print("Run=" + commons.PROG + " mode=" + mode + " alg=" + alg + " ds=" + ds + " threads=" + str(threads))
    for i in range(reps):
        stdout = subprocess.check_output(commons.PROG + " -a " + alg + " -d " + ds + " -n " + str(threads), cwd="../../",
                                         shell=True).decode('utf-8')
        print(str(i+1) + '/' + str(reps) + ': ' +
              mode + "," + alg + "," + ds + "," + str(threads))
        
        # Parse the logger output.
        for line in stdout.splitlines():
            # find throughput for the first graph
            # Use the "Throughput" in the enclave. Logger [ ENCL], [ DEBUG].
            if "Throughput" in line:
                throughput = re.findall("\d+\.\d+", line)[1]
                throughput_array.append(float(throughput))
            # find phase for the second graph
            # (phase meaning: "Phase Total", "Phase Build", "Phase Probe", etc.)
            if "Phase" in line:
                words = line.split()
                phase_name = words[words.index("Phase") + 1]
                value = int(re.findall(r'\d+', line)[-2])
                print (phase_name + " = " + str(value))
                if phase_name in dic_phases:
                    dic_phases[phase_name].append(value)
                else:
                    dic_phases[phase_name] = [value]

        print('Throughput = ' + str(throughput) + ' M [rec/s]')

    # Write the results into the csv data
    throughput = statistics.mean(throughput_array)
    s = (mode + "," + alg + "," + ds + "," + str(threads) + "," + str(round(throughput, 2)))
    f_throughput.write(s + '\n')
    f_throughput.close()

# Read the output csv file and plot the results.
def plot_throughput():

    # Read the output csv file
    plot_filename = "../img/Figure-04-Throughput-of-join-algorithms"
    csvf = open(fname_throughput, mode='r')
    csvr = csv.DictReader(csvf)
    all_data = list(csvr)
    # algos = ['CHT', 'PHT', 'PSM', 'RHT', 'RHO', 'RSM', 'NL', 'INL']
    algos = (set(map(lambda x:x['alg'], all_data)))
    datasets = sorted(set(map(lambda x:x['ds'], all_data)), reverse=True)
    modes = sorted(set(map(lambda x:x['mode'], all_data)))
    width = 0.4
    to_modes = [[y for y in all_data if y['mode'] == x] for x in modes]

    # graph per dataset
    plt.rc('axes', axisbelow=True)
    plt.rcParams.update({'font.size': 13})
    fig = plt.figure(figsize=(8,2.3))
    # plt.clf()
    legend_elements = [Patch(label='Hatches:', alpha=0),
                       Patch(facecolor='white', edgecolor='black',
                             hatch='\\\\', label='plain CPU'),
                       Patch(facecolor='white', edgecolor='black',
                             label='TEE')]

    to_datasets = [[y for y in all_data if y['ds'] == x] for x in datasets]
    for d in range(0, len(datasets)):
        ax = plt.subplot(1, 2, d+1)
        ax.yaxis.grid(linestyle='dashed')
        to_modes = [[y for y in to_datasets[d] if y['mode'] == x] for x in modes]
        for m in range(0, len(modes)):
            if m == 0:
                br = np.arange(len(algos))
                br = [x - 0.2 for x in br]
            else:
                br = [x + width for x in br]

            label = modes[m]
            # Because small place -> Maybe better use legends or short implementation names.
            algos = list(OrderedDict.fromkeys(list(map(lambda x:x['alg'],to_modes[m]))))
            colors = list(map(lambda x: commons.color_alg(x), algos)) # Change the algorithm colors in commons.
            hatch = '\\' if modes[m] == 'native' else ''
            plt.bar(br, list(map(lambda x: float(x['throughput']), to_modes[m])),
                    width=width, label=label, hatch=hatch, color=colors, edgecolor='black')
            for x, y in zip(br, list(map(lambda x: float(x['throughput']), to_modes[m]))):
                if y < 4:
                    plt.text(x-0.15, y+10, str(y), rotation=90, size=12)
            if d == 0:
                plt.ylabel("Throughput [M rec/s]", size=13)
            #plt.ylim([0, 230])
            plt.yticks(size=11)
            plt.xticks(np.arange(len(to_modes[m])), list(map(lambda x:x['alg'],to_modes[m])),
                       rotation=37, size=7)
            plt.title('(' + chr(97+d) + ") Dataset $\it{" + datasets[d] + "}$",
                      y=-0.4, size=13)

    fig.legend(handles=legend_elements, ncol=3, frameon=False,
               bbox_to_anchor=(0.048,0.83,1,0), loc="lower left",
               prop={'size': 12})
    commons.savefig(plot_filename + ".png", tight_layout=False)
# Main function
if __name__ == '__main__':

    #Load config file.
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, 'r:',['reps='])
    except getopt.GetoptError:
        print('Unknown argument')
        sys.exit(1)

    # config['reps'], config['experiment'] are variables from config file.
    for opt, arg in opts:
        if opt in ('-r', '--reps'):
            config['reps'] = int(arg)

    if config['experiment']:
        commons.remove_file(fname_throughput)
        commons.init_file(fname_throughput, "mode,alg,ds,threads,throughput\n")
        reps = 3 # config['reps']
        for mode in config['modes']:
            # Compile and build the enclave.
            commons.compile_app(mode, enclave_config_file='Enclave/Enclave2GB.config.xml')
            for ds in commons.get_test_dataset_names():
                for alg in ["CHT", "CHTopt", "RHT", "PRHO"]:
                    # Run the joins.
                    run_join(mode, alg, ds, config['threads'], reps ) 

    plot_throughput()
