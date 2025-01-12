# ScaWL

This repository contains supplementary code for the paper "ScaWL: Scaling $k$-WL (Weisfeiler-Leman) Algorithms in Memory and Performance on Shared and Distributed-Memory Systems".

You will need a version of MPI and OpenMP to compile this code. The set of instructions below are for compiling the code for maximum reproducibility and therefore are using
the same libraries that we used for running the test. The version of OpenMP used in this configuration is OpenMP 4.5.

The code was built with the Intel version of MPI (Intel(R) MPI Library 2019 Update 8 for Linux). GCC 7.3.0 was used as the compiler and Jemalloc 5.3.0 was used as the memory allocator. Assuming that both of these libraries are available as modules on your cluster follow the steps below:

1) Load the modules:

```bash
module load intel/20.2
module load jemalloc/5.3.0
```
2) Navigate to either the 2WL or 3WL directory, which ever version of ScaWL you want to use.

3) Build the scawl executable using the following command:

```bash
mpicxx -std=c++11 -I. -O3 -g -DEXCEL_OUTPUT -DALL_TO_ALL_V scawl.cpp -lpthread -ljemalloc -fopenmp -o scawl.exe
```
<br>

4) To prepare the graph data, navigate to the matrices folder and execute the script:
```bash
sh get_matrices.sh
```
This script will retrieve the matrices from the Florida University repository and perform all data preparation.
If the graphs cannot be connected to because of ssh certificate issues, you can temporarily bypass certificate validation
with the command:

```bash
sh get_matrices.sh no_cert_validation
```

This script will perform all the data preparation, including generating symmetric matrices when required as well as generating the
relabeled isomorphic matrices.


5) Once the executable is built and the matrix data has been retrieved, you can run all the jobs by navigating back to the executable folder
and executing the script:

```bash
sh run_all_jobs.sh
```

The batch jobs will be generated and an output file with the result of the isomorphism test will be produced for each test configuration.
The test configurations and graph set that is run correspond to the 2-WL and 3-WL tests in the research paper.
