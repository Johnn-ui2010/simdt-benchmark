#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# R:S tuples = 1:4

#import commons
import commons
import re
import statistics
import subprocess
import csv

filename = "../data/NL-output.csv"
mb_of_data = 131072


def run_join(prog, alg, size_r, size_s, threads, reps, mode):
    f = open(filename, "a")
    t_results = []
    s_results = []
    for i in range(0,reps):
        stdout = subprocess.check_output(prog + " -a " + alg + " -r " + str(size_r) + " -s " + str(size_s) + " -n " + str(threads), cwd="../../",shell=True).decode('utf-8')
        #print("mode "+ mode)

        for line in stdout.splitlines():
            if (mode == "native") and  ("Total join runtime" in line):
                time = re.findall("\d+\.\d+", line)[1]
                t = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(time))
                t_results.append(float(time))
                print (t)
            if (mode == "sgx") and ("throughput" in line):
                #print("I'm checking throughput.")
                throughput = re.findall("\d+\.\d+", line)[1]
                s = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(throughput))
                s_results.append(float(throughput))
                print (s)
            
    if(mode == "sgx"):
        time_res = 0
        # A googd idea is to change throughput to 4 decimal points in App.cpp.
        throughput_res = statistics.mean(s_results)
    else:
        throughput_res = 0
        time_res = statistics.mean(t_results)
    
    s = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(round(time_res,2)) + "," + str(throughput_res))
    print("AVG: " + s) 
    f.write(s + "\n")
    f.close()

if __name__ == '__main__':
    timer = commons.start_timer()
    total_sizes = [int(0.2*mb_of_data),   # 0.2 MB
               int(5 * mb_of_data), # 5 MB
               10 * mb_of_data,     # 10 MB = 1 310 720 tuples 
               20 * mb_of_data]      # 20 MB 
    reps = 3
    threads = 4
    modes = ["native","sgx"] # Place it later back
    commons.remove_file(filename)
    commons.init_file(filename, "mode,alg,threads,sizeR,sizeS,total_join_time,throughput\n")

    for mode in modes:
        #Compile the selected mode.
        commons.compile_app(mode)

        for tot_size in total_sizes:
            for alg in ["NL","NL_keys","NL_tuples"]: #commons.get_all_algorithms():
                r_size = tot_size * (1/5)
                s_size = tot_size * (4/5)
                run_join(commons.PROG, alg, r_size, s_size, threads, reps, mode)

    commons.stop_timer(timer)
