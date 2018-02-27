#!/bin/bash
#SBATCH -p serial
#SBATCH --mem=30GB
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
    
NAME="iii8"     # raw data folder name
PRESET="ndi305-40" # connection matrix and preset data file
FORMAT="png"        # figure format for saving
individOnly="false" # plot individual neurons’ figures only
statsOnly="true"   # plot populational figures only
ld="false"       # default to false, if already acquired sim data as ".mat" file, set to true
fit="false"     # fit the tuning curves for width.
rt=10        #   run time of simulation (exclude the warm-up time), some plots will give wrong firing rate if not set correctly as in adg_INPUT
npool=1           # number of cpus used to fit tuning curves
thres=0.25          # a gauge of firing rate threshold for active, non-active, non-active neurons were not included in some of the plots
NTHETA="12"         # number of angles in simulation
CONTRASTS="4"       # number of constrast used in simulation
MeanNorm="false"     # normalize before taking average

#matlab -nosplash -nodesktop -r "dbstop if error;plotFFT('$PRESET','$NAME',4,0,'png'); exit;"
#matlab -nosplash -nodesktop -r "plotISI('$NAME','$PRESET',4,5,$npool,'$FORMAT',$NTHETA); plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt,$fit,$npool);exit"
#matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt,$fit,$npool);exit"
#matlab -nosplash -nodesktop -r "CondVsAxis('$PRESET','$NAME',$NTHETA,$ld,'$FORMAT','$NAME');exit;"
NAME0='iii-nie'
#matlab -nosplash -nodesktop -r "plotFeedBackInh('$NAME0','$NAME','$PRESET',$NTHETA,'$FORMAT');exit;"
matlab -nosplash -nodesktop -r "plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt,$fit,$npool,$MeanNorm);exit"

echo Folder = $NAME
