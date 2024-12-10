# Benchmarking the Vectorized Join Implementations
This is a modified version of [TEEBench](https://github.com/agora-ecosystem/tee-bench) that was implemented by Maliszewski et al. We used this framework to evaluate our vectorized join implementations in SGX. The original paper is included [here](https://github.com/agora-ecosystem/tee-bench/blob/master/paper/What_Is_the_Price_for_Joining_Securely_Benchmarking_Equi-Joins_in_Trusted_Execution_Environments.pdf).

## Changes from original TEEBench 
We mostly did not change and retained the TEEBench
structure. General changes are listed here:
- Enabled the SIMD functionality.
- Added SIMD-vectorized implementations of CHT and NL join based on baseline implementations in TEEBench.
- Inserted PRHO implementation from [Balkesen](https://archive-systems.ethz.ch/node/334) et al. 
- Modified some original scripts for our experiment and added them to "simd_scripts" folder.
- Removed other unrelated original scripts.
- Removed the original paper, but retained the TEEBench instructions.
- Added diagrams and tables as results.
- Removed some unrelated joins (grace and stitch join).
- Edited the README file.

## Requirements
Make sure that your machine supports:
* Intel SGX v2.11 
* At least Ubuntu 18.04 LTS (experimented with Ubuntu 20.04.2 LTS)

### Installation
*  First, install the following dependencies:
```
  $ sudo apt update
  $ sudo apt install make gcc g++ libssl-dev python3-pip git-lfs  
  $ pip3 install matplotlib numpy pyyaml
```  
* Update SGX driver and install Intel SGX-SDK (and Microsoft SEAL).   
To simplify the process, we recommend using the script provided by TEEBench. 
```
$ cd simd_scripts 
$ ./install.sh
``` 
Troubleshooting: If your Intel SGX-SDK version does not work, install an older version.
* Now clone the repository using `git`:
```
$ git clone https://github.com/Johnn-ui2010/simdt-benchmark
``` 
* Next, download the IMDb dataset (from TEEBench). It can be downloaded from this tubCloud [repository](https://tubcloud.tu-berlin.de/s/KKXjkkS6GwD5pwn).
Add it to `simd_join_benchmarks/data` folder.

The chapter "Installation" in [TEE: How to Reproduce the Experimental Results](simd_scripts/TEEBench_How_to_Reproduce_the_Experimental_Results.pdf) describes the important installation steps again and can be used as reference. 

### SIMD libraries (with SSE and AVX2)
They are stored in `simd` folder.
It is not possible to just use "#include <immintrin.h>" in Intel SGX.  In SGX, the SIMD libraries should be included separately with "#include <mmintrin.h>" instead of "#include "mmintrin.h".

Thereafter, go to `simd_scripts` folder and check whether the installation was successful with this shell script. This script execute the vectorized CHT implementation `CHTfsimd` in SGX.
```
$ cd simd_scripts 
$ ./compile.sh
``` 
   
## Compilation
You can manually compile join algorithms in two modes: plain CPU or Intel SGX.

1. Intel SGX
 
   1. Run this command:  
   ` $ make clean && make -B sgx `
   2. Execute with:  
   ` $ ./app `
   
2. Plain CPU  
   
    1. Run this command:  
    ` $ make clean && make -B native CFLAGS=-DNATIVE_COMPILATION `  
    2. Execute with:  
    ` $ ./app `

## Execution
### Command line arguments
The following command line arguments are important for the experiment:
* `-a` - join algorithm name. Use: `CHT`, `CHTfsimd`, `RHT`, `PRHO`, `NL`, `NL_keys` or `NL_tuples`. Default: `RHT`
* `-d` - name of pre-defined dataset. Currently working: `cache-fit`, `cache-exceed`, `L`. Default: `none`
* `-n` - number of threads used to execute the join algorithm. Default: `2`
* `-r` - number of tuples of R relation. Default: `2097152`
* `-s` - number of tuples of S relation. Default: `2097152`
* `-t | --r-path` - filepath to build R relation. Default: `none`
* `-u | --s-path` - filepath to build S relation. Default `none`
* `-x` - size of R in MBs. Default: `none`
* `-y` - size of S in MBs. Default: `none`

One example is:
```
$ ./app -a CHTfsimd -r 1000000 -s 4000000 -n 4
```
## Scripts overview
The prefered way to execute the experiments is by using scripts.

We took some [TEEBench scripts](https://github.com/agora-ecosystem/tee-bench/tree/master/scripts) as baselines and modified them for our benchmark. In each script was commented, on which original script it was based on.
Our scripts are stored in `simd_scripts` folder.

Choose a script from `helper` folder and run it from command line, e.g.:
```
$ python3 exp1-off-the-shelf-performance_exp.py
```

Each script has an experiment only (exp) and plot only (plot) mode.
Replace `[A]` with exp or plot. 
Some experiments are just the same, one measuring the throughput and the other the join phases. 
We offer scripts that do it in one common experiment.
e.g. `exp1and9-off-the-shelf-experiment_exp.py`.

These experiments were considered in the evaluation chapter.

- `NL evaluation`
  - exp2-nl-scale-tuples-experiment_[A].py
  
- `Off-the-shelf-performance`: cache-fit vs cache-exceed data
(from TEEBench)`
  - exp1-off-the-shelf-performance_[A].py

  - exp9-phases-experiment_[A].py

- `Scaling the tuples for CHT by tables size` (R:S=1:4, 50 MB to 350 MB) 
  - exp3-cht-scale-tuples-experiment_[A].py

- `Real data experiment` (from TEEBench)
  - exp5-real-dataset-experiment_[A].py
  - exp10-real-dataset-experiment-phases_[A].py

- `Scaling the S-tuples for CHT by cache type` 
  (from TEEBench)
   - exp7-scale-s-epc-experiment_[A].py

- `Correctness test`
  - exp14-test-correctness_[alg].py , for alg = [NL, CHT, RHT]

## Evaluation
The results are stored in `data` and `img` folder. As mentioned before, each experiment has a plot scripts to visualize the results.
Additionally, you can look into this [Jupyter Notebook](/simd_join_benchmarks/simd_scripts/helpers/jup_notebook_tests.ipynb).
