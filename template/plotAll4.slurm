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

######## find variable description in pA.slurm or sbatch_sim4x3 ###########
echo NAME=$NAME
echo PRESET=$PRESET
echo FORMAT=$FORMAT
echo individOnly=$individOnly
echo statsOnly=$statsOnly
echo ld=$ld
echo fit=$fit
echo rt=$rt
echo npool=$npool
echo thres=$thres
echo fit=$fit
echo NTHETA=$NTHETA
echo CONTRASTS=$CONTRASTS

date
./assa

if [ -d "processed" ]; then
    rm -r 0125/ 0250/ 0500/ 1000/
    cp logNormalProfile.mat EPSPsize.mat
    
    mv processed/* ./
    rm -r processed
    
    cd ../
    
    #matlab -nosplash -nodesktop -r "dbstop if error;plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
    
    matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA,$individOnly,'$NAME',{},$thres,$statsOnly,false,$ld,$rt,$fit,$npool);exit"
    
    #matlab -nosplash -nodesktop -r "plotindividual('$NAME','$PRESET',[],'$FORMAT',$CONTRASTS,$NTHETA);exit"
    #matlab -nosplash -nodesktop -r "plotresult('$NAME','$PRESET',$NTHETA,$CONTRASTS,'$FORMAT');exit"
else
    echo failed
fi
echo Folder = $NAME
date
