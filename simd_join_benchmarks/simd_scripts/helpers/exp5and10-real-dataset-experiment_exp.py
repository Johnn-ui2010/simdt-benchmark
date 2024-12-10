#Modified the original script exp8-phases-experiment.py
#!/usr/bin/python3
import getopt
import subprocess
import re
import sys

import csv
import commons
import statistics
import yaml

fname_output = "../data/real-dataset-output.csv"
fname_phases = "../data/real-dataset-phases.csv"
r_path = 'data/imdb-name-basics.tsv'
s_path = 'data/imdb-title-principals.tsv'

def run_join(mode, alg,threads, reps):

    f = open(fname_output, 'a')
    f_phases = open(fname_phases, 'a')

    throughput_array = []
    throughput = ''
    dic_phases = {}
    print("Run=" + commons.PROG + " mode=" + mode + " alg=" + alg + " threads=" + str(threads))
    for i in range(reps):
        stdout = subprocess.check_output(commons.PROG + " -a " + alg + " -n " + str(threads)
                                         + " --r-path " + r_path + " --s-path " + s_path, cwd="../../",
                                         shell=True).decode('utf-8')
        print(str(i+1) + '/' + str(reps) + ': ' +
              mode + "," + alg + "," + str(threads))
        for line in stdout.splitlines():
            # find throughput for the first graph
            if "Throughput" in line:
                throughput = re.findall("\d+\.\d+", line)[1]
                throughput_array.append(float(throughput))

            # find phases for the second graph
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

        print('Throughput = ' + str(throughput) + ' M [rec/s]')

    #File for throughput
    throughput = statistics.mean(throughput_array)
    s = (mode + "," + alg + "," + str(threads) + "," + str(round(throughput, 2)))
    f.write(s + '\n')
    f.close()

    #File for the individual phases
    for x in dic_phases:
        res = statistics.mean(dic_phases[x])
        s = mode + "," + alg + "," + x + "," + str(res)
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

    if config['experiment']:
        commons.remove_file(fname_phases)
        commons.init_file(fname_phases, "mode,alg,phase,cycles\n")

        commons.remove_file(fname_output)
        commons.init_file(fname_output, "mode,alg,threads,throughput\n")

        for mode in config['modes']:
            commons.compile_app(mode)
            for alg in ["CHT","CHTfsimd","RHT","PRHO"]:
                run_join(mode, alg, config['threads'], config['reps'])

    commons.stop_timer(timer)