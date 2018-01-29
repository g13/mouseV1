#!/bin/bash
#SBATCH --mem-per-cpu=4GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --job-name=tmpPlot
#SBATCH -o tmp_%J.out

module load NYUAD/3.0
module load matlab
module load gcc/4.9.3
module list

ulimit -s unlimited

theme="isample40-t61l"
format='png';
loadData=true
tcReady=true

date

matlab -nosplash -nodesktop -r "simAssisVeff('$theme',$loadData,$tcReady,'$format','$theme');exit;"

date
