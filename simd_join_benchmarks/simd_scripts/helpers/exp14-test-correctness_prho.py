#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# Testing correctness of "PRHO_avx"
import commons
import getopt
import statistics

import subprocess
import re
import sys

import csv
from collections import OrderedDict
import yaml

filename = "../data/PRHO-correctness.csv"
r_path = 'data/imdb-name-basics.tsv'
s_path = 'data/imdb-title-principals.tsv'
mb_of_data = 131072


def run_join_normal(mode, prog, test_nr, alg, r_tuples, s_tuples, threads):
    f = open(filename, "a")
    "test_nr,alg,r_tuples,s_tuples,matches,is_correct,percent_diff\n"
    matches = 0

    #Real dataset
    if r_tuples == -1:
        stdout = subprocess.check_output(commons.PROG + " -a " + alg + " -n " + str(threads)
                                         + " --r-path " + r_path + " --s-path " + s_path, cwd="../../",
                                         shell=True).decode('utf-8')
    else:
        stdout = subprocess.check_output(prog + " -a " + alg + " -r " + str(r_tuples) + " -s " + str(s_tuples) + " -n " + str(threads), 
                                         cwd="../../",shell=True).decode('utf-8')
    
    for line in stdout.splitlines():
        if ("Matches" in line):
            matches = int(re.findall(r'\d+', line)[3])

        if ("Result tuples" in line) and (mode == "sgx"):
            matches = int(re.findall(r'\d+', line)[3])
            if matches == matches_correct:
                is_correct = True
            else:
                ratio_diff = (matches-matches_correct) / matches_correct

    s = (mode + "," + str(test_nr)  + "," + alg + "," + str(r_tuples) + "," + str(s_tuples) + "," + str(matches) + "," + "True" + "," + str(0) )
    f.write(s + "\n")
    f.close()

    return matches

def run_join_vect(mode,prog, test_nr, alg, r_tuples, s_tuples, threads,  matches_correct):
    f = open(filename, "a")

    matches = 0
    is_correct = False
    ratio_diff = 0

    #Real dataset
    if r_tuples == -1:
        stdout = subprocess.check_output(commons.PROG + " -a " + alg + " -n " + str(threads)
                                         + " --r-path " + r_path + " --s-path " + s_path, cwd="../../",
                                         shell=True).decode('utf-8')
    else:
        stdout = subprocess.check_output(prog + " -a " + alg + " -r " + str(r_tuples) + " -s " + str(s_tuples) + " -n " + str(threads), 
                                         cwd="../../",shell=True).decode('utf-8')
        
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
    total_sizes = [(1000000,1000000) ,  #Normal
                (405580,404564) ,     #Normal
               (10155,10031) , # Two numbers of tuples, not divisible by 4.
               (100440,100103) , # One number of tuples, not divisible by 4.
               (100103,100440) , # One number of tuples, not divisible by 4.

               (3042424,100000),  # R tuples bigger
               (1000000,3042424),  # S tuples bigger

                (10*mb_of_data,40*mb_of_data),#cache-fit
                (40*mb_of_data,10*mb_of_data), 
                (10*mb_of_data,120*mb_of_data), #bigger data set
                (1000,400000), #bigger tuple difference
                (10012,4000650) #bigger tuple difference
               ]      
    reps = 1
    threads = 2
    commons.remove_file(filename)
    commons.init_file(filename, "mode,test_nr,alg,r_tuples,s_tuples,matches,is_correct,ratio_diff\n")
    normal_alg = "RHT"
    #Compile the selected mode.
    for mode in config['modes']:
    #for mode in ['native']:
        commons.compile_app(mode)
        i = 0
        for tot_size in total_sizes:
            r_tuples = tot_size[0]
            s_tuples = tot_size[1]
            
            matches_correct = run_join_normal(mode,commons.PROG, i,normal_alg, r_tuples, s_tuples, threads)

            for alg in ["PRHO","PRHO_avx"]: #commons.get_all_algorithms():
                run_join_vect(mode,commons.PROG,i, alg, r_tuples, s_tuples, threads, matches_correct)

            i=i+1
    
        #Real dataset
        matches_correct = run_join_normal(mode,commons.PROG, i,normal_alg, -1, -1, threads)
        for alg in ["PRHO","PRHO_avx"]: #commons.get_all_algorithms():
            run_join_vect(mode,commons.PROG,i, alg, -1, -1, threads, matches_correct)

    print("Please go to csv file for results.")
    commons.stop_timer(timer)



