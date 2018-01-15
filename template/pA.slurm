#!/bin/bash
#SBATCH -p serial
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --time=0:35:00
#SBATCH -o pA_%A.out
#SBATCH --job-name=plotAll
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END
module purge
module load matlab
module load gcc/4.9.1
module list
ulimit -s unlimited
    
NAME="sample40-6t"     # raw data folder
PRESET="1xu-ndi305-40-s911.mat" # connection matrix and preset data file
CONTRASTS="4"       # # of constrast used in simulation
FORMAT="png"        # figure format for saving
NTHETA="12"         # number of angles in simulation
individOnly="false" # plot individual neurons’ figures only
thres=0.45          # a gauge of firing rate threshold for active, non-active
statsOnly="false"    # plot populational figures only
rt=5        #   run time of simulation
ld="false"  #   default to false, if already acquired sim data as ".mat" file, set to true

matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt);exit"
#matlab -nosplash -nodesktop -r "plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly);exit"
#matlab -nosplash -nodesktop -r "startFit(10800,'$NAME',1:$CONTRASTS,12,12,true);ppp;exit"
#matlab -nosplash -nodesktop -r "dbstop if error;plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
#matlab -nosplash -nodesktop -r "dbstop if error;plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
#matlab -nosplash -nodesktop -r "plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
#matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');exit"

#matlab -nosplash -nodesktop -r "startFit(10800,'$NAME',1:$CONTRASTS,12,12,true);plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,12);exit"
#matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly);exit"

#matlab -nosplash -nodesktop -r "dbstop if error;plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
#matlab -nosplash -nodesktop -r "dbstop if error;plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
#matlab -nosplash -nodesktop -r "plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"

echo Folder = $NAME
