#!/bin/sh
#SBATCH --constraint=haswell,e5-2660-v3-2.60ghz
#SBATCH --output=%x.%j.out                # Output file (stdout)
#SBATCH --error=%x.%j.err                 # Error file (stderr)

# Get the parameters passed by sbatch
directory1=$DIRECTORY1
file1=$FILE1
directory2=$DIRECTORY2
file2=$FILE2
num_cores=$CORES
echo "The value of directory1 is $directory1"
echo "The value of file1 is $file1"
echo "The value of directory2 is $directory2"
echo "The value of file2 is $file2"
echo "The value of cores is $num_cores"
echo $PATH
echo Running test..
export MLX5_SINGLE_THREADED=0
export I_MPI_JOB_TIMEOUT=432000
export LD_PRELOAD=/opt/apps/jemalloc/5.3.0/lib/libjemalloc.so
export GOMP_CPU_AFFINITY="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores
export OMP_STACKSIZE=1G
unset I_MPI_PMI_LIBRARY
export OMP_NUM_THREADS=$num_cores
mpirun --bind-to none ./scawl.exe $directory1/$file1 $directory2/$file2
