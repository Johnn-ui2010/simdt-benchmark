import re
import statistics
import subprocess
import csv
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

filename = "data/simd_file2.csv"

def plot(operation):

    print("The plots are created...")

    out_file = ["images/"+ operation + "_native.pdf","images/"+ operation + "_sgx.pdf",
            "images/"+ operation + "_native_throughput.pdf","images/"+ operation + "_sgx_throughput.pdf",
            "images/"+ operation + "_speedup.pdf"]
    
    size_names = [
            '${size}$ < L2',
            'L2 < ${size}$ < L3',
            'L3 < ${size}$ < EPC',
            'EPC < ${size}$'
        ]
    colors_list = {"scalar":"#ff5454", "unroll":"#ff9900", "simd_sse":"#0e8052", "simd_avx":"#a13905",
                   "simd_sse_native":"#c2ffba", "simd_avx_native":"#b89967", "simd_sse_sgx":"#0e8052","simd_avx_sgx":"#a13905",
                   }
    
    df = pd.read_csv(filename)

    df_native = df.loc[df['mode'] == 'native']
    df_sgx = df.loc[df['mode'] == 'sgx']

    df_native_op = df_native.loc[df['operation'] == operation]
    df_sgx_op = df_sgx.loc[df['operation'] == operation]

    #Calculate normal execution times in [ms]
    y_native_op = {}
    for t in df_native_op['type'].unique():
        y_native_op[t] = ( df_native_op.loc[df_native['type'] == t]["time"].reset_index(drop=True).to_list() )
    y_native_op = pd.DataFrame(y_native_op)

    y_sgx_op = {}
    for t in df_sgx_op['type'].unique():
        y_sgx_op[t] = ( df_sgx_op.loc[df_sgx['type'] == t]["time"].reset_index(drop=True).to_list() )
    y_sgx_op = pd.DataFrame(y_sgx_op)

    #Calculate execution time per tuple times.
    mb_of_data = 262144
    num_tuples = [int(0.2*mb_of_data),int(6.4*mb_of_data),16.0*mb_of_data,100.0*mb_of_data]
    
    y_native_op_throughput = y_native_op.copy()

    for i in range(len(num_tuples)):  
        y_native_op_throughput.loc[i,:] = 1000 * num_tuples[i] / 1000000 / y_native_op_throughput.loc[i,:] 

    y_sgx_op_throughput = y_sgx_op.copy()

    for i in range(len(num_tuples)):  
        y_sgx_op_throughput.loc[i,:] = 1000 * num_tuples[i] / 1000000 / y_sgx_op_throughput.loc[i,:] 

    fig, ax = plt.subplots(1,1)

    # algorithm categories/implementations
    type_cat = df_native_op['type'].unique()

    x =np.arange( 0,len(size_names) )

    #y_speedup = {"unroll_native":[], "simd_sse_native":[],"simd_avx_native":[],"unroll_sgx":[],"simd_sse_sgx":[],"simd_avx_sgx":[]}
    y_speedup = {"scalar":[1.0, 1.0, 1.0, 1.0], "simd_sse_native":[],"simd_avx_native":[],"simd_sse_sgx":[],"simd_avx_sgx":[]}

    i = 0
    for y in [y_native_op, y_sgx_op]:
        ax = y.plot(kind="bar", color=colors_list)
        ax.legend(type_cat)
        ax.set_xticks(x, size_names, rotation = 0)
        ax.set_xlabel("Total table size",fontsize=14)
        ax.set_ylabel("CPU execution time in ms",fontsize=14)
        offset = -0.2
        for t in df_native_op['type'].unique():
            if i == 0:
                for m in range(4):
                    if y[t][m] < 10:
                        ax.text(m + offset, y[t][m]+1, str(round(y[t][m],2)), rotation=90, size=8)
                    
                    if (t == "simd_sse") or (t == "simd_avx"):
                        y_speedup[t + "_native"].append( round(y["scalar"][m]/y[t][m],2) )
                    #if t == "vectorized":
                    #   txt = str(round(y["scalar"][m]/y["vectorized"][m],2))
                    #   ax.text(m + offset, y[t][m]+2, txt, rotation=0, size=10,color="red")
                    
                offset += 0.15

            else:
                for m in range(4):
                    if y[t][m] < 2:
                        ax.text(m + offset, y[t][m]+0.7, str(round(y[t][m],2)), rotation=90, size=8)
                    if (t == "simd_sse") or (t == "simd_avx"):
                        y_speedup[t + "_sgx"].append( round(y["scalar"][m]/y[t][m],2) )
                    #if t == "vectorized":
                    #    txt = str(round(y["scalar"][m]/y["vectorized"][m],2))
                    #    ax.text(m + offset, y[t][m]+2, txt, rotation=0, size=10,color="red")
                    
                offset += 0.15

        plt.savefig(out_file[i], transparent = False, bbox_inches = 'tight', dpi=200)
        print("The plots is successfully created and is in ", out_file[i])
        i = i+1
        ax.clear()
    
    for y in [y_native_op_throughput,y_sgx_op_throughput]:
        ax = y.plot(kind="bar",color=colors_list)
        ax.legend(type_cat)
        ax.set_xticks(x, size_names, rotation = 0)
        ax.set_xlabel("Total table size",fontsize=12)
        ax.set_ylabel("Throughput in M rec/s",fontsize=12)
        plt.savefig(out_file[i], transparent = False, bbox_inches = 'tight', dpi=200)
        print("The plots is successfully created and is in ", out_file[i])
        ax.clear()
        i = i+1



    #Calculate and plot SIMD speedup 
    print(y_speedup)
    y_speedup = pd.DataFrame(y_speedup)
    if operation == 'ADD':
        y_speedup.plot(color=colors_list, marker='x',figsize=(6, 4))
        plt.legend(fontsize = 8, loc='upper right')
    else:
        y_speedup.plot(color=colors_list, marker='x',legend=False,figsize=(6, 4))

    

    plt.xticks(x, size_names, rotation = 0,fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Total table size",fontsize=12)
    plt.ylabel("SIMD speedup",fontsize=12)
    plt.ylim([0, 3.6])
    plt.savefig(out_file[4], transparent = False, bbox_inches = 'tight')#, dpi=200)
    print("The plots is successfully created and is in ", out_file[4])
    
if __name__ == '__main__':
    #Enter the operation.
    operations = ['ADD','AND','SHIFT_RIGHT']
    for op in operations:
        plot(op)