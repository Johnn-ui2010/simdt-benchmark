# Microbenchmarking the operators
This is the code for our microbenchmark.
The three operators (ADD, AND and SHIFT_RIGHT) are part of our experiment.

We wrote four versions for each operator. The first version is the scalar version. In the second, we unrolled the loop. In version 3, the operators were vectorized with the SSE instruction set and in the fourth with AVX.

We modified the `SampleEnclave` example from the initial Intel SGX repository to fit our experiment.

## Installation
Please follow the steps from our [other benchmark](/../simd_teebench/simd_join_benchmarks/README.md) for SGX.

### SIMD libraries (with SSE and AVX2)
They are stored in "simd" folder.
It is not possible to just use "#include <immintrin.h>" in Intel SGX.  In SGX, the SIMD libraries should be included separately with "#include <mmintrin.h>" instead of "#include "mmintrin.h".

Thereafter, go to `MathsExperiments` folder and check whether the installation was successful with this shell script.
```
$ cd scripts
$ ./compile_sgx.sh
```

## Experiment
### Running the experiment 
We assume that you are in `Microbenchmark` folder.
#### Running it manually
Plese run the experiments in SGX with the following command.

```
$ make SGX_PRERELEASE=1 SGX_DEBUG=0
```

For plain mode use:
```
$ make native
```

#### Using scripts
Please go to `scripts` folder and execute the following Python script.
```
$ cd scripts
$ python3 script_exp.py
```

## Evaluation
The results are available as a csv file in `data` folder and as diagrams in `img` folder.
The diagrams can be plotted using a Python script, too.
```
$ cd scripts
$ python3 script_plot.py
```
Alternatively, we provided a [Jupyter Notebook](/simd_microbenchmarks/MathsExperiments/scripts/notebook.ipynb).