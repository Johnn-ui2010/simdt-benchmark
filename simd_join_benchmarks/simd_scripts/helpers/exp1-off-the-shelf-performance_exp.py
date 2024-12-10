#!/usr/bin/python3
#  Original script modified: exp1-off-the-shelf-performance.py

import getopt
import subprocess
import re
import sys

import csv
import commons
import statistics
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

# Main function
if __name__ == '__main__':
    timer = commons.start_timer()

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
        reps = config['reps']
        for mode in config['modes']:
        #for mode in ["sgx"]:
            # Compile and build the enclave.
            commons.compile_app(mode)
            for ds in commons.get_test_dataset_names():
                for alg in ["CHT", "CHTfsimd", "RHT", "PRHO"]:
                    # Run the joins.
                    run_join(mode, alg, ds, config['threads'], reps ) 
    commons.stop_timer(timer)