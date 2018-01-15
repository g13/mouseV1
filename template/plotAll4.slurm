#!/bin/bash
#SBATCH -p serial 
#SBATCH --mem=40GB
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=12
#SBATCH -J pA4
#SBATCH -o pA4_%J.out
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END

module purge
module load matlab
module load gcc/4.9.1
ulimit -s unlimited

#NAME="4-M-21"
NAME="sample40-6t"
#PRESET="1xu-nd35-35-s911.mat"
PRESET="1xu-ndi305-40-s911.mat"
#PRESET="1xu-nmLDoO-of32-test-s911.mat"
CONTRASTS="4"
FORMAT="png"
NTHETA="12"
thres=0.1
individOnly="false"
statsOnly="false"
fit="true"
rt=5
ld="false"
fit="false"

date
./assa

if [ -d "processed" ]; then
    cp generateConMat* gg-* profile.out conMat.mat processed
    cp logNormalProfile.m plotAll4.slurm sbatch_sim4x3 cMatrix_summary adg_INPUT sample_id.dat processed
    cp logNormalProfile.mat processed/EPSPsize.mat
    cp s4*.out processed
    cp pA4*.out processed
    
    mv V1.* processed/
    mv processed ../
    
    rm -r *
    mv ../processed/* ./
    cd ../
    rm -r processed
    
    #matlab -nosplash -nodesktop -r "startFit(10800,'$NAME',1:$CONTRASTS,12,12,true);exit"
    #matlab -nosplash -nodesktop -r "dbstop if error;plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
    
    if [ "$fit" == 'true' ]; then
        matlab -nosplash -nodesktop -r "startFit(10800,'$NAME',1:$CONTRASTS,12,12,true);plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt);exit"
    else
        matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt);exit"
    fi
    
    ##matlab -nosplash -nodesktop -r "plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
    ##matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');exit"
    #
else
    echo failed
fi
echo Folder = $NAME
date
