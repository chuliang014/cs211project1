#!/bin/bash
#SBATCH -o mydata.out
#SBATCH -e mydata.err
#SBATCH -J mydata
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=11000
#SBATCH --time=6:0:0
# make sure you are using gcc-4.7.2
module load gcc/gcc-4.7.2

# compile your code on target node
gcc -o main main.c -lm

# &> is overwrite update
./main &> lb1_data.txt

# >> is add-on update
