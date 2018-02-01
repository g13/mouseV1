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
module load gcc/4.9.3
ulimit -s unlimited

############# Please find variable description in pA.slurm ###############
NAME="sample40-neo6"
PRESET="1xu-ndi305-40-s911.mat"
FORMAT="png"
individOnly="false"
statsOnly="true"
ld="false"
fit="true"
rt=10
npool=12
thres=0.25
fit="true"
NTHETA="12"
CONTRASTS="4"

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
    
    matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt,$fit,$npool);exit"
    
    ##matlab -nosplash -nodesktop -r "plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
    ##matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');exit"
    #
else
    echo failed
fi
echo Folder = $NAME
date
