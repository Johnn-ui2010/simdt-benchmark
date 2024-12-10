#!/usr/bin/python3

# Original script modified: exp7-scale-s-experiment-SUDO.py
# Similar experiment like in Barber CHT paper.
# Fix the number of R tuples and vary the S tuples. And then the same in opposite.

import getopt
import subprocess
import re
import sys

import csv
import commons
import statistics
import yaml

filename_fixed_r_phases = "../data/scale_expon_fixed_r_phases.csv"
filename_fixed_s_phases = "../data/scale_expon_fixed_s_phases.csv"
mb_of_data = 131072


def run_join(prog, alg, size_r, size_s, fixed_r, threads, reps, mode):
    if (fixed_r):
        f_phases = open(filename_fixed_r_phases, 'a')
    else:
        f_phases = open(filename_fixed_s_phases, 'a')

    s_results = []
    throughput = ''
    dic_phases = {}
    
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

            # find phases for the second graph (phases)
            if "Phase" in line:
                words = line.split()
                phase_name = words[words.index("Phase") + 1]
                value = int(re.findall(r'\d+', line)[-2])
                print (phase_name + " = " + str(value))
                if phase_name in dic_phases:
                    dic_phases[phase_name].append(value)
                else:
                    dic_phases[phase_name] = [value]
            if "Build (cycles)" in line:
                #print(line)
                value = int(re.findall(r'\d+', line)[3])
                print("Build" + " = " + str(value))
                if "Build" in dic_phases:
                    dic_phases["Build"].append(value)
                else:
                    dic_phases["Build"] = [value]
            if "Probe (cycles)" in line:
                #print(line)
                value = int(re.findall(r'\d+', line)[3])
                print("Probe" + " = " + str(value))
                if "Probe" in dic_phases:
                    dic_phases["Probe"].append(value)
                else:
                    dic_phases["Probe"] = [value]

    throughput_res = statistics.mean(s_results)
    s = (mode + "," + alg + "," + str(threads) + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + str(round(throughput_res,2)))
    print("AVG: " + s) 

    #File for the individual phases
    for x in dic_phases:
        res = statistics.mean(dic_phases[x])
        s = mode + "," + alg + "," + str(round(size_r/mb_of_data,2)) + "," + str(round(size_s/mb_of_data,2)) + "," + x + "," + str(res)
        f_phases.write(s + '\n')

    f_phases.close()

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
    for opt, arg in opts:
        if opt in ('-r', '--reps'):
            config['reps'] = int(arg)

    reps = config['reps'] 
    threads = config['threads']
    modes = ["native"]
    
    commons.remove_file(filename_fixed_r_phases)
    commons.remove_file(filename_fixed_s_phases)
    commons.init_file(filename_fixed_r_phases, "mode,alg,threads,sizeR,sizeS,phase,cycles\n")
    commons.init_file(filename_fixed_s_phases, "mode,alg,threads,sizeR,sizeS,phase,cycles\n")
    #
    for mode in modes:

       #Compile the selected mode.
        commons.compile_app(mode)

        for i in range(0,10):
            # More R
            for alg in ["CHT", "CHTfsimd", "RHT", "PRHO"]:
                r_size = mb_of_data * 2 ** (10) 
                s_size = mb_of_data * 2 ** (i)

                run_join(commons.PROG, alg, r_size, s_size, 1, threads, reps, mode)

            # More S
            for alg in ["CHT", "CHTfsimd", "RHT", "PRHO"]:
                r_size = mb_of_data * 2 ** (i)
                s_size = mb_of_data *2 ** (10)
                run_join(commons.PROG, alg, r_size, s_size, 0, threads, reps, mode)

    commons.stop_timer(timer)
