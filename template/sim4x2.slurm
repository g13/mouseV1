#!/bin/bash
#SBATCH -p serial
#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH -J sim4x2
#SBATCH -o s42_%J.out



module purge
module load gcc/4.9.3


export OMP_NUM_THREADS=12
export OMP_STACKSIZE=5G
ulimit -s unlimited

date
./ompaSep2 4
date
