#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --job-name=gg
#SBATCH -o gg_%J.out
#SBATCH --mail-user=wd554@nyu.edu
#SBATCH --mail-type=END
module purge
module load matlab
fdr='test'              # target simulation folder
theme='ndi305-40'       # simulation name
lgnfile='1xu-'$theme'-s911' # data file name for lgn position and map to v1
TYPE='l06'              # additional type name
conmat=conMat-$theme'_'$TYPE  # name for connection matrix
coMa='coMa-60x60-'$lgnfile   # name for RF pearson coefficient matrix
EPSPprofile='logNormal'     # type of distribution for EPSP
edist='coGauss'             # E-E connection strength type
eiSpecific='true'          # E<-I connection specificity over orientation preference
antiphase='false'       # antiphase in terms of on-off subregion, so that RF coefficient is bigger if on-off cross match between two neurons
reciprocal=0.0          # constraint on percentage of reciprocal connections 0.0 for no constraint
drawInh='true'      # if draw inhibitory example; 
ld='false'          # re-generate lgn to v1 map
reverse=2           # different types weight contribution see RFcoeff_O_beta.m

tw=0.25             # weight of contributions from orientation preference and RFcoeff to connection preference, dependent on "reverse", examine in RFcoeff_O_beta.m
temper=30            # randomness to temper with the preset orientation of recruiting ellipses
ntheta=12           # number of angles simulated
ex=10               #  # of example plots
p2s=40.0            # estimated converting scale from connection strength to resulting PSP
mu=0.45             # average EPSP in mV
ieSpecific='false'  # keep it false, its coMat is not implemented in RFcoeff_O_beta.m (set I<-E specific connection over oriention preference)
sig=1.16           # standard deviation for a lognormal distribution
mu0=0.45           # with a mean value of mu0
nsig=3             # cutoff EPSP histogram at # of std 
nbins=30            # number of bins for EPSP distribution
spread='false'      # unrelated test variable, assign EPSP value to E->E connections by subregion normalized distance, ignore.
lowestRatio=0.02    # related to variable "spread", ignore.
if [ "$ld" == 'false' ]; then
   matlab -nosplash -nodesktop -r "gMap2;exit"
fi
# pre copy
cp generateConMat.m $fdr/generateConMat-$theme'_'$TYPE.m
cp gg $fdr/gg-$theme'_'$TYPE
cp $lgnfile.mat $fdr/
cp lgnmap.out-$lgnfile $fdr/lgnmap.out-$lgnfile
cp logNormalProfile.m $fdr/logNormalProfile-$theme'_'$TYPE.m
cp lgn2v1map/$lgnfile/gMap2.m $fdr/gMap2-$theme'_'$TYPE.m
cp lgn2v1map/$lgnfile/lgn2v1Map_beta.m $fdr/lgn2v1Map_beta-$theme'_'$TYPE.m

date

#cp RFcoeff_V.m $fdr/RFcoeff_O_beta-$theme'_'$TYPE.m
#matlab -nosplash -nodesktop -r "dbstop if error;RFcoeff_V('$lgnfile',60,60,3,true,true,$tw,true,'png',$ld,$temper,$reverse);logNormalProfile(400,rand(400,1),$nbins,$mu,$sig,true);generateConMat('$lgnfile','t2','gauss',$reciprocal,'coGauss',911,true,'png','n','logNormal','$conmat',$eiSpecific,'$coMa');checkConnection_EI('mee-$conmat','$lgnfile','$coMa',randi(8640,1,$ex),9,'png',true,$ntheta,$drawInh);checkConnection_EI('mee-$conmat','$lgnfile','$coMa',1:8640,8,'png',true,$ntheta,$drawInh);exit"
#cp RFcoeff_O.m $fdr/RFcoeff_O-$theme'_'$TYPE.m
cp RFcoeff_O_beta.m $fdr/RFcoeff_O_beta-$theme'_'$TYPE.m
matlab -nosplash -nodesktop -r "dbstop if error;logP.nbins=$nbins;logP.mu=$mu;logP.sigma=$sig;logP.mu0=$mu0;logP.nsig=$nsig;logP.lR=$lowestRatio;logP.p2s=$p2s;logP.spread=$spread;RFcoeff_O_beta('$lgnfile',60,60,3,true,true,$tw,true,'png',$ld,$temper,$reverse,false,$ieSpecific,$antiphase);generateConMat('$lgnfile','t2','gauss',$reciprocal,'$edist',911,true,'png','n','$EPSPprofile','$conmat',$eiSpecific,$ieSpecific,'$coMa',logP);checkConnection_EI('mee-$conmat','$lgnfile','$coMa',randi(8640,1,$ex),9,'png',true,$ntheta,$drawInh,$spread,'$TYPE');checkConnection_EI('mee-$conmat','$lgnfile','$coMa',1:8640,8,'png',true,$ntheta,$drawInh,$spread,'$TYPE');exit"

cp $coMa'_more.mat' $fdr/coMat-$theme'_'$TYPE.mat
cp $conmat.mat $fdr

cp logNormalProfile.mat $fdr/logNormalProfile-$theme'_'$TYPE.mat
cp logNormalProfile.out $fdr/profile-$theme'_'$TYPE.out
date