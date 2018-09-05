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
    matlab -nosplash -nodesktop -r "checkConnection_EI('mee-$conmat','$theme','$coMa',1:8640,8,'png',true,$ntheta,$drawInh,$spread,'$TYPE');exit"
    
    #cp $coMa'_more.mat' ../$fdr/coMat-$theme'_'$TYPE.mat
    #cp $conmat.mat ../$fdr
    #
    #cp logNormalProfile.mat ../$fdr/logNormalProfile-$theme'_'$TYPE.mat
    #cp logNormalProfile.out ../$fdr/profile-$theme'_'$TYPE.out
date
