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

colors_list = {"NL":"#ff5454", "NL_keys":"#0e8052", "NL_tuples":"#a13905",
            "NL_keys_native":"#c2ffba", "NL_tuples_native":"#b89967", 
            "NL_keys_sgx":"#0e8052","NL_tuples_sgx":"#a13905",}

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

    ax = y_native.plot(kind = "bar",color = colors_list)
    ax.legend(alg_cat)
    ax.set_xticks(x, x_tuples_size, rotation = 0)
    ax.set_xlabel("Total tables size in MB",)
    ax.set_ylabel("Time in s")
    ax.set_title("native mode",fontsize = 11)
    print(y_native["NL"][0])
    print(x_tuples_size[0])

    # See the original script: exp1-off-the-shelf-performance.py
    offset = -0.2
    y_speedup = {"NL":[1.0, 1.0, 1.0, 1.0],"NL_tuples_native":[], "NL_tuples_sgx":[]}

    for alg in csvf_native['alg'].unique():
        for i in range(4):
            if y_native[alg][i] < 10:
                ax.text(i + offset, y_native[alg][i]+10, str(y_native[alg][i]), rotation=90, size=9)
            if alg == "NL_tuples":
                y_speedup["NL_tuples_native"].append( round(y_native["NL"][i]/y_native["NL_tuples"][i],3) )

                txt = str(round(y_native["NL"][i]/y_native["NL_tuples"][i],2))
                ax.text(i + offset, y_native[alg][i]+50, txt, rotation=0, size=12,color="red")
                    
        offset += 0.15

    commons.savefig('../img/Fig-e2-NL-output_native.pdf')

    y_native_throughput = y_native.copy()
    mb_of_data = 131072
    sizes = [0.2*mb_of_data,5*mb_of_data,10*mb_of_data,20*mb_of_data]
    for i in range(len(sizes)):  
        y_native_throughput.loc[i,:] = sizes[i] / 1000000 / y_native_throughput.loc[i,:] 

    x =np.arange( 0,len(x_tuples_size) )

    
    ax = y_native_throughput.plot(kind="bar", color = colors_list)
    ax.legend(alg_cat)
    ax.set_xticks(x, x_tuples_size, rotation = 0)
    ax.set_ylim([0,3])
    ax.set_xlabel("Total tables size in MB",)
    ax.set_ylabel("Throughput in M rec/s")
    plt.title("native mode",fontsize = 11)

    #L2, L3 cache lines
    ax.vlines(x=0.5, ymin= 0, ymax= 1.5, color = "blue")
    ax.text(x=0.4, y= 1.6, s = "<-L2", color = "blue",fontsize=11)
    ax.vlines(x=2.5, ymin= 0, ymax= 1.5, color = "blue")
    ax.text(x=2.4, y= 1.6, s = "<-L3", color = "blue",fontsize=11)

    offset = -0.2
    for alg in csvf_native['alg'].unique():
        for i in range(4):
            if y_native_throughput[alg][i] < 0.25:
                ax.text(i + offset, y_native_throughput[alg][i]+0.05, str(round(y_native_throughput[alg][i],3)), rotation=90, size=9)
            if alg == "NL_tuples":
                if i == 0:
                    txt = str(round(y_native_throughput["NL_tuples"][i]/y_native_throughput["NL"][i],2))
                    ax.text(i + offset, y_native_throughput[alg][i]+0.2, txt, rotation=0, size=12,color="red")
                else:
                    txt = str(round(y_native_throughput["NL_tuples"][i]/y_native_throughput["NL"][i],2))
                    ax.text(i + offset, y_native_throughput[alg][i]+0.5, txt, rotation=0, size=12,color="red")
        offset += 0.15
    commons.savefig('../img/Fig-e2-NL-output_throughput_native.pdf')

    ax.clear()
    #Plot the SGX mode
    fig, ax = plt.subplots(1,1)
    alg_cat = csvf_sgx['alg'].unique()
    x =np.arange( 0,4 )

    ax = y_sgx.plot(kind = "bar", color = colors_list)
    ax.legend(alg_cat)
    ax.set_xticks(x, x_tuples_size, rotation = 0)
    ax.set_xlabel("Total tables size in MB")
    ax.set_ylabel("Throughput in M rec/s")
    ax.set_ylim([0,3])
    ax.set_title("SGX mode",fontsize = 11)

    #L2, L3 cache lines
    ax.vlines(x=0.5, ymin= 0, ymax= 1.5, color = "blue")
    ax.text(x=0.4, y= 1.6, s = "<-L2", color = "blue",fontsize=11)
    ax.vlines(x=2.5, ymin= 0, ymax= 1.5, color = "blue")
    ax.text(x=2.4, y= 1.6, s = "<-L3", color = "blue",fontsize=11)
    
    # See the original script: exp1-off-the-shelf-performance.py
    offset = -0.2
    for alg in csvf_sgx['alg'].unique():
        for i in range(4):
            if y_sgx[alg][i] < 0.25:
                ax.text(i + offset, y_sgx[alg][i]+0.05, str(round(y_sgx[alg][i],3)), rotation=90, size=9)
            if alg == "NL_tuples":
                y_speedup["NL_tuples_sgx"].append( round(y_sgx["NL_tuples"][i]/y_sgx["NL"][i],3) )

                if i == 0:
                    txt = str(round(y_sgx["NL_tuples"][i]/y_sgx["NL"][i],2))
                    ax.text(i + offset, y_sgx[alg][i]+0.2, txt, rotation=0, size=12,color="red")
                else:
                    txt = str(round(y_sgx["NL_tuples"][i]/y_sgx["NL"][i],2))
                    ax.text(i + offset, y_sgx[alg][i]+0.5, txt, rotation=0, size=12,color="red")

        offset += 0.15
    commons.savefig('../img/Fig-e2-NL-output_sgx.pdf')

    print(y_speedup)
    y_speedup = pd.DataFrame(y_speedup)
    #{'simd_sse_native': [1.9, 1.55, 1.54, 1.52], 'simd_avx_native': [3.52, 2.79, 2.71, 2.61], 'simd_sse_sgx': [1.04, 0.76, 0.81, 0.87], 'simd_avx_sgx': [2.48, 1.32, 1.34, 1.32]}
    y_speedup.plot(color=colors_list)
    plt.xticks(x, x_tuples_size, rotation = 0)
    plt.xlabel("Total tables size in MB")
    plt.ylabel("SIMD speedup")

    commons.savefig('../img/Fig-e2-NL-output_speedup.pdf')
    
if __name__ == '__main__':
    plot()
