#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# filename_more_r: R:S tuples = 1:4
# filename_more_s: R:S tuples = 4:1

#[TO DO, not correct.]
#import commons
import commons
import re
import statistics
import subprocess
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filename_more_r = "../data/big_scale_more_r.csv"
filename_more_s = "../data/big_scale_more_s.csv"

out_file = ["../img/Fig-e8-big-scale_native_more_r_line.pdf",'../img/Fig-e8-big-scale_native_more_r_bar.pdf',
                "../img/Fig-e8-big-scale_native_more_s_line.pdf",'../img/Fig-e8-big-scale_native_more_s_bar.pdf',
]

graph_titles = ["native mode (R:S = 4:1)","native mode (R:S = 4:1)",
                "native mode (R:S = 1:4)","native mode (R:S = 1:4)"]
mb_of_data = 131072

def plot():
    i = 0
    for f in [filename_more_r,filename_more_s]:
        csvf = pd.read_csv(f)
        csvf = csvf.loc[csvf['alg'] != 'PRHO']
        csvf_native = csvf.loc[csvf['mode'] == 'native']

        csvf_native["total_size"] = csvf_native["sizeR"] + csvf_native["sizeS"]
        x_tuples_size = csvf_native.loc[csvf_native['alg'] == 'CHT']["sizeR"] + csvf_native.loc[csvf_native['alg'] == 'CHT']["sizeS"]
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
        x =np.arange( 0,len(x_tuples_size) )


        ax = y_native.plot()
        ax.legend(alg_cat)
        ax.set_xticks(x, x_tuples_size, rotation = 0)
        ax.set_xlabel("Total tables size in MB",)
        ax.set_ylabel("Throughput in M rec/s")
        ax.set_title(graph_titles[i],fontsize = 11)
        commons.savefig(out_file[i])
        i= i+1
        ax.clear()

        ax = y_native.plot(kind="bar")
        ax.legend(alg_cat,prop={'size': 6})
        ax.set_xticks(x, x_tuples_size, rotation = 0)
        ax.set_xlabel("Total tables size in MB",)
        ax.set_ylabel("Throughput in M rec/s")
        ax.set_title(graph_titles[i],fontsize = 11)
        commons.savefig(out_file[i])
        i= i+1
        ax.clear()

if __name__ == '__main__':
    plot()