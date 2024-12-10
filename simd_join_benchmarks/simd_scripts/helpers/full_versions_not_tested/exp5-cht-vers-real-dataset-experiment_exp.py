#!/usr/bin/python3
# Modified the original: exp5-real-dataset-experiment.py
import getopt
import subprocess
import re
import sys

import csv
import commons
import statistics

import yaml


fname_output = "../data/real-dataset-output.csv"
r_path = 'data/imdb-name-basics.tsv'
s_path = 'data/imdb-title-principals.tsv'


def run_join(mode, alg, threads, reps):

    f = open(fname_output, 'a')

    throughput_array = []
    throughput = ''
    print("Run=" + commons.PROG + " mode=" + mode + " alg=" + alg + " threads=" + str(threads))
    for i in range(reps):
        stdout = subprocess.check_output(commons.PROG + " -a " + alg + " -n " + str(threads)
                                         + " --r-path " + r_path + " --s-path " + s_path, cwd="../../",
                                         shell=True).decode('utf-8')
        print(str(i+1) + '/' + str(reps) + ': ' +
              mode + "," + alg + "," + str(threads))
        for line in stdout.splitlines():
            # find throughput
            if "Throughput" in line:
                throughput = re.findall("\d+\.\d+", line)[1]
                throughput_array.append(float(throughput))

        print('Throughput = ' + str(throughput) + ' M [rec/s]')

    throughput = statistics.mean(throughput_array)
    s = (mode + "," + alg + "," + str(threads) + "," + str(round(throughput, 2)))
    f.write(s + '\n')
    f.close()

if __name__ == '__main__':
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, 'r:',['reps='])
    except getopt.GetoptError:
        print('Unknown argument')
        sys.exit(1)
    for opt, arg in opts:
        if opt in ('-r', '--reps'):
            config['reps'] = int(arg)

    if config['experiment']:
        commons.remove_file(fname_output)
        commons.init_file(fname_output, "mode,alg,threads,throughput\n")

        for mode in config['modes']:
            #commons.compile_app(mode, enclave_config_file='Enclave/Enclave2GB.config.xml', debug=True)
            #commons.compile_app(mode, enclave_config_file='Enclave/Enclave2GB.config.xml')
            commons.compile_app(mode)
            for alg in ["CHT","CHTfsimd"]:
                run_join(mode, alg, config['threads'], config['reps'])

