# ScaWL: Scaling k-WL (Weisfeiler-Leman) Algorithms in Memory and Performance on Shared and Distributed-Memory Systems

This repository contains supplementary code for the paper "ScaWL: Scaling $k$-WL (Weisfeiler-Leman) Algorithms in Memory and Performance on Shared and Distributed-Memory Systems".

You will need a version of MPI and OpenMP to compile this code. The set of instructions below is for compiling the code for maximum reproducibility and therefore is using
the same libraries that we used for running the test. The version of OpenMP used in this configuration is OpenMP 4.5.

The code was built with the Intel version of MPI (Intel(R) MPI Library 2019 Update 8 for Linux). GCC 7.3.0 was used as the compiler and Jemalloc 5.3.0 was used as the memory allocator. Assuming that both of these libraries are available as modules on your cluster follow the steps below:

1) Load the modules:

```bash
module load intel/20.2
module load jemalloc/5.3.0
```

2) Navigate to either the 2WL or 3WL directory, whichever version of ScaWL you want to use.

3) Navigate to the file run_tests.sh. This is the lowest-level script file that launches the ScaWL executable and contains the environment variables that the executable runs with. Make sure the LD_PRELOAD environment variable is pointing to the appropriate jemalloc library. By default, it is pointing to the location on the cluster ScaWL was tested on. Inside run_tests.sh, it would look something like this:

```bash
export LD_PRELOAD=/path/to/your/jemalloc/5.3.0/lib/libjemalloc.so
```

4) Build the ScaWL executable using the following command:

```bash
mpicxx -std=c++11 -I. -O3 -DEXCEL_OUTPUT -DALL_TO_ALL_V scawl.cpp -lpthread -ljemalloc -fopenmp -o scawl.exe
```

5) To prepare the graph data, navigate to the matrices folder and execute the script:
```bash
sh get_matrices.sh
```
This script will retrieve the matrices from the University of Florida repository and perform all data preparation.
If the graphs cannot be connected to because of SSH certificate issues, you can temporarily bypass certificate validation
with the command:

```bash
sh get_matrices.sh no_cert_validation
```

This script will perform all the data preparation, including generating symmetric matrices when required as well as generating the
relabeled isomorphic matrices.


6) Once the executable is built and the matrix data has been retrieved, you can run all the jobs by navigating back to the executable folder
and executing the script:

```bash
sh run_all_jobs.sh
```

The batch jobs will be generated and an output file with the result of the isomorphism test will be produced for each test configuration.
The test configurations and graph set that is run correspond to the 2-WL and 3-WL tests in the research paper. The output files will have the format
\<matrixName\>_\<nodeCount\>_\<coreCount\>.\<jobNumber\>.out. For example: 'appu.mtx_1Nodes_1Core.55898693.out'. Inside the result file, there will be an excel formatted
result in the format \<graph1Location\>,\<graph2Location\>,\<numVertices\>,\<worldSize\>,\<omp_get_num_procs()\>,\<omp_get_max_threads()\>,\<timeInSeconds\>,\<isIsomorphic\>.
The scripts were developed to make sure that the requested number of cores and the numbers of threads (by setting OMP_NUM_THREADS=cores) are the same. Therefore, in your result
file, these two numbers should match. The example excel result from the output file above is 'appu/appu.mtx,appu/appu.mtx,14000,1,1,1,4236.13,true'. Below is table defining the Excel output in more detail:

| Term        | Definition                                                                     |
|-------------|--------------------------------------------------------------------------------|
| **graph1Location**  | The directory location and name of graph 1.                            |
| **graph2Location**  | The directory location and name of graph 2.                            |
| **numVertices**  | The number of vertices contained in both graph 1 and graph 2.             |
| **worldSize**  | The MPI world size. In ScaWL, this is equal to the number of compute nodes. |
| **omp_get_num_procs()**  | The number of cores per compute node.                             |
| **omp_get_max_threads()**  | The number of threads per compute node.                         |
| **timeInSeconds**  | The time it takes to finish the isomorphism test in seconds.            |
| **isIsomorphic**  | The result of the isomorphism test.                                      |


## Building without Jemalloc 5.3.0

If Jemalloc 5.3.0 is not available on your cluster, try contacting your cluster administrator to get it installed. If that is not possible, the steps to build without Jemalloc are as follows:

1) Load just the Intel module:

```bash
module load intel/20.2
```

2) Navigate to either the 2WL or 3WL directory, whichever version of ScaWL you want to use.

3) Navigate to the file run_tests.sh. Remove the line with the LD_PRELOAD environment variable entirely.

4) Build the ScaWL executable without the command to link to Jemalloc:

```bash
mpicxx -std=c++11 -I. -O3 -DEXCEL_OUTPUT -DALL_TO_ALL_V scawl.cpp -lpthread -fopenmp -o scawl.exe
```

You will now be using the standard memory allocator. Since the default allocator is less performant than Jemalloc, you will likely get slower overall speeds and more memory fragmentation. This fragmentation will likely lead to more frequent failed allocations. Therefore, it is highly recommended to use Jemalloc if your goal is to reproduce the results
of the research paper as closely as possible. However, the standard allocator will get you up and running.


# Appendix

This appendix covers some of the details of the scripts and programs that are used to both generate the batch jobs on the cluster as well as prepare the graph data.

## Data Preparation Scripts

| Code File        | Definition                                                                     |
|-------------|--------------------------------------------------------------------------------|
| **degree_calc.c**  | This program was used by the research paper to calculate min and max degree for all graphs.                            |
| **bidirection_adder.c**  | This program generates the other half of symmetric graphs in the matrix market format.                            |
| **iso_create.c**  | This program creates an isomorphic version of the input graph by reversing the labels.             |
| **get_matrices.h**  | This script retrieves all the matrices for the University of Florida website, compiles bidirection_adder.c and iso_create.c, and finally calls make_iso.h or make_sym_and_iso.h. Whether make_iso.h or make_sym_and_iso.h is called depends on the requirements of the graph.                            |
| **make_iso.h**  | This script calls the iso_create program to genenate and isomorphic version of the input graph. It appends _iso.mtx to the mtx file.         |
| **make_sym_and_iso.h**  | This scripts calls both the bidirection_adder and iso_create program back-to-back, appending both _sym_iso.mtx to the mtx file.             |


## Scripts that Launch ScaWL

The scripts in the 2-WL and 3-WL directory that launch the different configuations of ScaWL are arranged in a hiearchy that ultimately passes several parameters to the bottom-level script 'run_test.sh'. Below is a table describing what each of the different scripts do.

| Code File        | Definition                                                                     |
|-------------|--------------------------------------------------------------------------------|
| **run_all_graphs.sh**  | This script calls the top-most script for each graph, passing both graphs and the time to allot for the batch job.                            |
| **run_tests_entire_set.sh**  | This script calls the scripts for both multnode and single node test configurations, passing down graph and time allotment parameters for the batch job.                            |
| **run_tests_multnode_driver.sh**  | This script calls all the lower-level scripts to run all multinode configurations, passing down graph and time allotment parameters for the batch job, in addition to adding parameters for memory, node count, and core count.             |
| **run_tests_singlenode_driver.sh**  | This script calls all the lower-level scripts to run all single node configuratino, passing down graph and time allotment parameters for the batch job, in addition to adding parameters for memory, node count, and core count.  |
| **run_tests_main.sh**  | This script passes down all accumulated parameters and some environment variables to the final script run_test.sh.         |
| **run_test.sh**  | This script takes all the parameters that have been passed and uses them to launch ScaWL in the configuration the parameters specify.             |