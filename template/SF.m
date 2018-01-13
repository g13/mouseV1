fdr0 = {'mE4-12-'};
sfLabel = {'20','30','60'};
principalSF = 2;
contrast = 1:4;
ntheta = 4;
format = 'png';
outputfdr ='SF-pairs'; 
lgnfile = '1xu-n3-s911.mat';

fdr = strcat(fdr0,sfLabel);
neuronlist = -1; % get neuronlist and return
ranking = {};
%format = 'psc2';
nSF = length(fdr);

%for i = 1:nSF
%	plotresult(fdr{i},lgnfile,ntheta,contrast,format,outputfdr);
%end
%for i = 1:nSF
%	j = seq(i);
	[ranking, neuronlist] = plotind1eps(fdr{principalSF},lgnfile,neuronlist,format,max(contrast),ntheta,true,outputfdr,ranking);
%end

spatialFreqTC(fdr,sfLabel,lgnfile,neuronlist,ranking,format,contrast,outputfdr);
disp('finished');

%fdr = {'t2SF-1','t2SF-2','t2SF-4','t2SF-8','t2SF-16'};
%fdr = {'tnSF-4'};
%neuronlist = [];
%ranking = {};
%for i = 1:length(fdr)
%	plotresult(fdr{i},lgnfile,ntheta,contrast,format,fdr{1});
%end
%for i = 1:nSF
%	j = seq(i);
%	[ranking, neuronlist] = plotindividual(fdr{j},lgnfile,neuronlist,format,contrast,ntheta,true,outputfdr,ranking);
%end
%
%%spatialFreqTC(fdr,sfLabel,lgnfile,neuronlist,format,contrast,outputfdr);
%disp('finished t2SF');
%
%neuronlist = [];
%ranking = {};
%for i = 1:nSF
%	j = seq(i);
%	[ranking, neuronlist] = plotindividual(fdr{j},lgnfile,neuronlist,format,contrast,ntheta,true,outputfdr,ranking);
%end
%fdr = {'tsSF-4'};
%for i = 1:nSF
%	j = seq(i);
%	[ranking, neuronlist] = plotindividual(fdr{j},lgnfile,neuronlist,format,contrast,ntheta,true,outputfdr,ranking);
%end
