#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# R:S tuples = 1:4

import commons
import re
import statistics
import subprocess
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
            if (mode == "native") and ("Total join runtime" in line):
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
        throughput_res = statistics.mean(s_results)
    else:
        throughput_res = 0
        time_res = statistics.mean(t_results)
    
    s = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(round(time_res,2)) + "," + str(round(throughput_res,2)))
    print("AVG: " + s) 
    f.write(s + "\n")
    f.close()


def plot():
    #Read the CSV
    csvf = pd.read_csv(filename)
    csvf_native = csvf.loc[csvf['mode'] == 'native']
    csvf_sgx = csvf.loc[csvf['mode'] == 'sgx']
    csvf_native["total_size"] = csvf_native["sizeR"] + csvf_native["sizeS"]

    x_tuples_size = csvf_native.loc[csvf_native['alg'] == 'NL']["sizeR"] + csvf_native.loc[csvf_native['alg'] == 'NL']["sizeS"]
    x_tuples_size = x_tuples_size.reset_index(drop=True).astype('str')

    #Get the total runtime in native mode
    y_native = {}
    for alg in csvf_native['alg'].unique():
        y_native[alg] = ( csvf_native.loc[csvf_native['alg'] == alg]["total_join_time"].reset_index(drop=True).to_list() )
    y_native = pd.DataFrame(y_native)

    #Get the throughput in SGX mode
    y_sgx = {}
    for alg in csvf_sgx['alg'].unique():
        y_sgx[alg] = ( csvf_sgx.loc[csvf_sgx['alg'] == alg]["throughput"].reset_index(drop=True).to_list() )
    y_sgx = pd.DataFrame(y_sgx)

    #Plot the native mode
    fig, ax = plt.subplots(1,1)
    alg_cat = csvf_native['alg'].unique()
    x =np.arange( 0,4 )

    ax = y_native.plot(kind = "bar")
    ax.legend(alg_cat)
    ax.set_xticks(x, x_tuples_size, rotation = 0)
    ax.set_xlabel("Total tables size in MB",)
    ax.set_ylabel("Time in s")
    ax.set_title("native mode",fontsize = 11)

    commons.savefig('../img/NL-output_native.png')

    ax.clear()
    #Plot the SGX mode
    fig, ax = plt.subplots(1,1)
    alg_cat = csvf_sgx['alg'].unique()
    x =np.arange( 0,4 )

    ax = y_sgx.plot(kind = "bar")
    ax.legend(alg_cat)
    ax.set_xticks(x, x_tuples_size, rotation = 0)
    ax.set_xlabel("Total tables size in MB")
    ax.set_ylabel("Throughpt in M rec/s")
    ax.set_title("SGX mode",fontsize = 11)
    commons.savefig('../img/NL-output_sgx.png')


if __name__ == '__main__':
    timer = commons.start_timer()
    total_sizes = [int(1*mb_of_data),   # 1 MB
               int(5 * mb_of_data), # 1 MB
               10 * mb_of_data,     # 10 MB = 1 310 720 tuples 
               20 * mb_of_data]      # 20 MB 
    reps = 3
    threads = 2
    modes = ["native","sgx"]
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

    plot()
    commons.stop_timer(timer)
