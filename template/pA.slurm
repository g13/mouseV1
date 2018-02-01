#!/bin/bash
#SBATCH -p serial
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH -o pA_%A.out
#SBATCH --job-name=plotAll
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END

module purge
module load matlab
module load gcc/4.9.3
module list
ulimit -s unlimited
    
NAME="sample40-neo6-adap"     # raw data folder name
PRESET="1xu-ndi305-40-s911.mat" # connection matrix and preset data file
FORMAT="png"        # figure format for saving
individOnly="false" # plot individual neurons’ figures only
statsOnly="true"   # plot populational figures only
ld="false"  #   default to false, if already acquired sim data as ".mat" file, set to true
fit="true"  # fit the tuning curves for width.
rt=5        #   run time of simulation
npool=12
thres=0.5          # a gauge of firing rate threshold for active, non-active, non-active neurons were not included in some of the plots
NTHETA="12"         # number of angles in simulation
CONTRASTS="4"       # # of constrast used in simulation

matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt,$fit,$npool);exit"

echo Folder = $NAME
