#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# Similar experiment like in Barber CHT paper.
# Fix the number of R tuples and vary the S tuples. And then the same in opposite.

#import commons
import commons
import re
import statistics
import subprocess
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filename_fixed_r = "../data/scale_expon_fixed_r.csv"
filename_fixed_s = "../data/scale_expon_fixed_s.csv"

out_file = ["../img/Fig-e4-output_native_fixed_r.pdf","../img/Fig-e4-output_native_fixed_s.pdf"]

graph_titles = ["native mode: Fixed #R = 1024 MB = 10 ^(5), vary #S","native mode: Fixed #S = 1024 MB = 10 ^(5), vary #R"]

mb_of_data = 131072

def plot():
    i = 0
    for f in [filename_fixed_r,filename_fixed_s]:
        csvf = pd.read_csv(f)
        csvf = csvf.loc[csvf['alg'] != 'PRHO']
        csvf_native = csvf.loc[csvf['mode'] == 'native']

        if csvf_native.loc[csvf_native['alg'] == 'CHT']["sizeR"][0] > csvf_native.loc[csvf_native['alg'] == 'CHT']["sizeS"][0]:
            x_tuples_size = csvf_native.loc[csvf_native['alg'] == 'CHT']["sizeS"]
        else:
            x_tuples_size = csvf_native.loc[csvf_native['alg'] == 'CHT']["sizeR"]
        x_tuples_size = x_tuples_size.reset_index(drop=True)

        #Get the throughput in native mode
        y_native = {}
        # There are a lot of CHTversions. Please don't choose all for better readibility.
        for alg in csvf_native['alg'].unique():
            y_native[alg] = ( csvf_native.loc[csvf_native['alg'] == alg]["throughput"].reset_index(drop=True).to_list() )
        y_native = pd.DataFrame(y_native)

        fig, ax = plt.subplots(1,1)

        # algorithm categories/implementations
        alg_cat = csvf_native['alg'].unique()
        print(alg_cat)
        x =np.arange( 0,len(x_tuples_size) )


        ax = y_native.plot()
        ax.legend(alg_cat)
        ax.set_xticks(x, x_tuples_size, rotation = 90)
        ax.set_xlabel("Total tables size in MB",)
        ax.set_ylabel("Throughput in M rec/s")
        ax.set_title(graph_titles[i],fontsize = 11)
        commons.savefig(out_file[i])
        i= i+1
        ax.clear()

if __name__ == '__main__':
    plot()