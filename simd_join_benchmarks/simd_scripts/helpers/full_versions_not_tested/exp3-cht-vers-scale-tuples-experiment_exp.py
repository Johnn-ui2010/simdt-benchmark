#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# filename_more_r: R:S tuples = 1:4
# filename_more_s: R:S tuples = 4:1

import commons
import re
import statistics
import subprocess
import csv

filename_more_r = "../data/CHT_versions_more_r.csv"
filename_more_s = "../data/CHT_versions_more_s.csv"
mb_of_data = 131072


def run_join(prog, alg, size_r, size_s, more_r, threads, reps, mode):
    if (more_r):
        f = open(filename_more_r, "a")
    else:
        f = open(filename_more_s, "a")

    t_results = []
    s_results = []
    for i in range(0,reps):
        stdout = subprocess.check_output(prog + " -a " + alg + " -r " + str(size_r) + " -s " + str(size_s) + " -n " + str(threads), cwd="../../",shell=True).decode('utf-8')
        #print("mode "+ mode)

        for line in stdout.splitlines():
            if ("Throughput" in line):
                #print("I'm checking throughput.")
                throughput = re.findall("\d+\.\d+", line)[1]
                s = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(throughput))
                s_results.append(float(throughput))
                print (s)
            

    throughput_res = statistics.mean(s_results)

    
    s = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(round(throughput_res,2)))
    print("AVG: " + s) 
    f.write(s + "\n")
    f.close()

if __name__ == '__main__':
    timer = commons.start_timer()

    reps = 3
    threads = 2
    modes = ["native","sgx"]
    
    commons.remove_file(filename_more_r)
    commons.remove_file(filename_more_s)
    commons.init_file(filename_more_r, "mode,alg,threads,sizeR,sizeS,throughput\n")
    commons.init_file(filename_more_s, "mode,alg,threads,sizeR,sizeS,throughput\n")
    #
    for mode in modes:

       #Compile the selected mode.
        commons.compile_app(mode)

        for i in range(50,351,50):
            # More R
            for alg in ["CHT","CHTPa1","CHTPa3","CHTb1","CHTb1s", "CHTb2", "CHTb2s","CHTopt"]:
                tot_size = i * mb_of_data
                r_size = tot_size * (4/5)
                s_size = tot_size * (1/5)
                run_join(commons.PROG, alg, r_size, s_size, 1, threads, reps, mode)

            # More S
            for alg in ["CHT","CHTPr1","CHTPr2","CHTopt"]:
                tot_size = i * mb_of_data
                r_size = tot_size * (1/5)
                s_size = tot_size * (4/5)
                run_join(commons.PROG, alg, r_size, s_size, 0, threads, reps, mode)

            

    commons.stop_timer(timer)
