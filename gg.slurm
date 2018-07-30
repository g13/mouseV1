#!/bin/bash
#SBATCH --mem=30G
#SBATCH --time=0:40:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=gg
#SBATCH -o gg_%J.out
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END
module purge
module load matlab/R2015b
set -e
date
matlab -nosplash -nodesktop -r "gMap2;exit"
cp lgnmap.out-$theme ../$fdr/
cp $theme.mat ../$fdr/
date

if [ ! "$testLGNonly" == true ]; then
    matlab -nosplash -nodesktop -r "logP.nbins=$nbins;logP.mu=$mu;logP.sigma=$sig;logP.mu0=$mu0;logP.nsig=$nsig;logP.lR=$lowestRatio;logP.p2s=$p2s;logP.spread=$spread;RFcoeff_O_beta('$theme','$TYPE',60,60,3,true,$tw,true,'png',$ldRF,$temper,$temperGrad,$ieSpecific,$antiphase,$threads,$wC);generateConMat('$theme','t2','gauss',$reciprocal,'$edist',911,true,'png','n','$EPSPprofile','$conmat',$eiSpecific,$ieSpecific,'$coMa',logP,$eiSig,$eeSig,$ieSig);checkConnection_EI('mee-$conmat','$theme','$coMa',randi(8640,1,$ex),9,'png',true,$ntheta,$drawInh,$spread,'$TYPE');checkConnection_EI('mee-$conmat','$theme','$coMa',1:8640,8,'png',true,$ntheta,$drawInh,$spread,'$TYPE');exit"
    
    cp $coMa'_more.mat' ../$fdr/coMat-$theme'_'$TYPE.mat
    cp $conmat.mat ../$fdr
    
    cp logNormalProfile.mat ../$fdr/logNormalProfile-$theme'_'$TYPE.mat
    cp logNormalProfile.out ../$fdr/profile-$theme'_'$TYPE.out
    date
fi
