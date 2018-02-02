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

echo $NAME
echo $PRESET
echo $FORMAT
echo $individOnly
echo $statsOnly
echo $ld
echo $fit
echo $rt
echo $npool
echo $thres
echo $fit
echo $NTHETA
echo $CONTRASTS

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
