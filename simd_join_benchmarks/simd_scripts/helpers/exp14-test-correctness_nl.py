#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py

import commons
import getopt
import statistics

import subprocess
import re
import sys

import csv
from collections import OrderedDict
import yaml

filename = "../data/NL-correctness.csv"
mb_of_data = 131072


def run_join_normal(mode, prog, test_nr, alg, r_tuples, s_tuples, threads):
    f = open(filename, "a")
    "test_nr,alg,r_tuples,s_tuples,matches,is_correct,percent_diff\n"
    matches = 0

    stdout = subprocess.check_output(prog + " -a " + alg + " -r " + str(r_tuples) + " -s " + str(s_tuples) + " -n " + str(threads), cwd="../../",shell=True).decode('utf-8')
    for line in stdout.splitlines():
        if ("Matches" in line):
            matches = int(re.findall(r'\d+', line)[3])
        if ("Result tuples" in line) and (mode == "sgx"):
            matches = int(re.findall(r'\d+', line)[3])
    
    s = (mode + "," + str(test_nr)  + "," + alg + "," + str(r_tuples) + "," + str(s_tuples) + "," + str(matches) + "," + "True" + "," + str(0) )
    f.write(s + "\n")
    f.close()

    return matches

def run_join_vect(mode, prog, test_nr, alg, r_tuples, s_tuples, threads,  matches_correct):
    f = open(filename, "a")

    matches = 0
    is_correct = False
    ratio_diff = 0
    stdout = subprocess.check_output(prog + " -a " + alg + " -r " + str(r_tuples) + " -s " + str(s_tuples) + " -n " + str(threads), cwd="../../",shell=True).decode('utf-8')
    for line in stdout.splitlines():
        if ("Matches" in line):
            matches = int(re.findall(r'\d+', line)[3])
            if matches == matches_correct:
                is_correct = True
            else:
                ratio_diff = (matches-matches_correct) / matches_correct

        if ("Result tuples" in line) and (mode == "sgx"):
            matches = int(re.findall(r'\d+', line)[3])
            if matches == matches_correct:
                is_correct = True
            else:
                ratio_diff = (matches-matches_correct) / matches_correct

    s = (mode + "," + str(test_nr) + "," + alg + "," + str(r_tuples) + "," + str(s_tuples) + "," + str(matches) + "," + str(is_correct) + "," + str(ratio_diff) )
    f.write(s + "\n")
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

    # config['reps'], config['experiment'] are variables from config file.
    for opt, arg in opts:
        if opt in ('-r', '--reps'):
            config['reps'] = int(arg)

    timer = commons.start_timer()
    total_sizes = [(100000,100000) ,  #Normal
                (45580,44564) ,     #Normal
               (1155,1031) , # Two numbers of tuples, not divisible by 4.
               (1440,1103) , # One number of tuples, not divisible by 4.
               (1103,1440) , # One number of tuples, not divisible by 4.

               (342424,10000),  # R tuples bigger
               (100000,342424),  # S tuples bigger
              
               ]      
    reps = 1
    threads = 2
    commons.remove_file(filename)
    commons.init_file(filename, "mode,test_nr,alg,r_tuples,s_tuples,matches,is_correct,ratio_diff\n")
    normal_alg = "NL"

    for mode in config['modes']:
        #Compile the selected mode.
        commons.compile_app(mode)
        i = 0
        for tot_size in total_sizes:
            r_tuples = tot_size[0]
            s_tuples = tot_size[1]
            
            matches_correct = run_join_normal(mode, commons.PROG, i,normal_alg, r_tuples, s_tuples, threads)

            for alg in ["NL_keys","NL_tuples"]: #commons.get_all_algorithms():
                run_join_vect(mode, commons.PROG,i, alg, r_tuples, s_tuples, threads, matches_correct)

            i=i+1
    print("Please go to csv file for results.")
    commons.stop_timer(timer)


