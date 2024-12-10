#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# Original script modified: exp6-scale-s-experiment-SUDO.py
import commons
import getopt
import statistics

import subprocess
import re
import sys

import csv
from collections import OrderedDict
import yaml

filename = "../data/scale-r-output.csv"
mb_of_data = 131072

def run_join(prog, alg, size_r, size_s, threads, reps, mode):
    f = open(filename, "a")
    results = []
    for i in range(0,reps):
        stdout = subprocess.check_output(prog + " -a " + alg + " -r " + str(size_r) + " -s " + str(size_s) + " -n " + str(threads), cwd="../../",shell=True).decode('utf-8')

        for line in stdout.splitlines():
            if "Throughput" in line:
                throughput = re.findall("\d+\.\d+", line)[1]
                r = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(throughput))
                results.append(float(throughput))
                print (r)
    res = statistics.mean(results)
    s = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(round(res,2)))
    print ("AVG : " + r)
    f.write(s + "\n")
    f.close()


if __name__ == '__main__':
    timer = commons.start_timer()

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

    max_r_size_mb = 128
    s_sizes = [int(0.2*mb_of_data),   # ~205kB
               int(6.4 * mb_of_data), # 6.4 MB
               32 * mb_of_data,       # 20 MB
               100 * mb_of_data]      # 100 MB
    mode = "sgx"
    commons.compile_app(mode)
    commons.remove_file(filename)
    commons.init_file(filename, "mode,alg,threads,sizeR,sizeS,throughput\n")
    #
    for s_size in s_sizes:
        for alg in ["CHT","CHTfsimd"]:
            for i in range(8, max_r_size_mb+1, 32):
                run_join(commons.PROG, alg, i*mb_of_data, s_size, config["threads"], config["reps"], mode)

    commons.stop_timer(timer)
