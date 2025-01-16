#!/bin/sh
echo $PATH
sbatch -N $1 --time=$8 --job-name=3WL_$5_$1Nodes_$3Core --mem=$2 -p kamiak --cpus-per-task=$3 --export=DIRECTORY1=$4,FILE1=$5,DIRECTORY2=$6,FILE2=$7,CORES=$3,PATH=$PATH,LD_LIBRARY_PATH=$LD_LIBRARY_PATH run_test.sh
