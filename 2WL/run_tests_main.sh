echo $PATH
if [ -z "$6" ] || [ -z "$7" ]; then
    sbatch -N $1 --time=1-12:00:00 --job-name=$5_$1Nodes_$3Core --mem=$2 -p kamiak --cpus-per-task=$3 --export=DIRECTORY1=$4,FILE1=$5,DIRECTORY2=$4,FILE2=$5,CORES=$3,PATH=$PATH,LD_LIBRARY_PATH=$LD_LIBRARY_PATH run_tests.sh
else
    sbatch -N $1 --time=1-12:00:00 --job-name=$5_$1Nodes_$3Core --mem=$2 -p kamiak --cpus-per-task=$3 --export=DIRECTORY1=$4,FILE1=$5,DIRECTORY2=$6,FILE2=$7,CORES=$3,PATH=$PATH,LD_LIBRARY_PATH=$LD_LIBRARY_PATH run_tests.sh
fi