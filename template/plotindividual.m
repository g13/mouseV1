function [ranking, neuronlist] = plotindividual(theme,lgn,neuronlist,format,contrastLevel,ntheta,individOnly,outputfdr,ranking,thres,statsOnly,neuronlistOnly,dataRead,rt,fit,npool)
% position =[500,150,1000,650];set(0, 'OuterPosition', position);
% position = [50,50,1400,900];set(0, 'DefaultFigurePosition', position);
addpath(genpath('../matlab_Utilities'));
gL = 50;
vrest = 0;
vexcit = 35/11;
vinhib = -0.4;
disp(['gL = ',num2str(gL)]);
disp(['vE = ',num2str(vexcit)]);
disp(['vI = ',num2str(vinhib)]);
disp(['vR = ',num2str(vrest)]);
LineWidth = 2;
set(groot,'defaultLineLineWidth',LineWidth);
set(groot,'defaultErrorbarLineWidth',LineWidth);
FontSize = 14;
set(groot,'defaultAxesFontSize',FontSize);
set(groot,'defaultTextFontSize',FontSize);
LegendOffset = 1;
set(groot,'defaultLegendFontSize',FontSize-LegendOffset);
pPosition = [0, 0, 1280, 720];
if nargin < 16
    npool = 12;
    if nargin < 15
        fit = true;
        if nargin < 14
            rt = 5;
            if nargin < 13
                dataRead = false;
                if nargin < 12
                    neuronlistOnly = false;
                    if nargin < 11
                        statsOnly = false;
                        if nargin < 10
                            thres = 0.0; % gauge for active
                            if nargin < 9
                                ranking = {};
                                if nargin < 8
                                    outputfdr = theme;
                                    if nargin < 7
                                        individOnly = false;
                                        if nargin < 6
                                            ntheta = 12;      
                                            if nargin < 5
                                                contrastLevel = 4;
                                                if nargin < 4                  
                                                    format = '';
                                                    if nargin < 3
                                                        neuronlist = [];
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if ~isempty(format)
    if strcmp(format,'psc2')
        printDriver = ['-de',format];
        format = 'eps';
    else
        printDriver = ['-d',format];
    end
    dpi = '-r100';
end
meanTimeLine = true;
%meanTimeLine = false;
% current = true;
current = false;
%operiod = true;
operiod = false;
%dThetaSize = true;
dThetaSize = false;
%Phase = true;
Phase = false;

ndperiod = 25;
polarplot = 0;
pOSI = true;
rng('shuffle');

    load(lgn);
    nLGN = nLGN(:,1);
    ntotal = p.nv1;
    CVsortedID = zeros(ntotal,contrastLevel);
    PKsortedID = CVsortedID;
    conLabel = cell(contrastLevel,1);
    load([theme,'/EPSPsize.mat']);
if ~dataRead
    tic;
    for i=1:contrastLevel
        disp(['loading c',num2str(i)]);
        eps = 125*2^(i-1);
        conLabel(i) = {[num2str(eps/10,'%.1f'),'%']};
    %     eps = 250*i;
        DIR = [theme,'/',num2str(eps,'%04d')];
            [nP(i),tC(i),aC] = readDataAll(DIR,ntheta,ntotal,ndperiod);
        pC(i).Vstd=zeros(ntotal,ndperiod);
        pC(i).gtotstd=zeros(ntotal,ndperiod);
        pC(i).Itotstd=zeros(ntotal,ndperiod);
        pC(i).Veffstd=zeros(ntotal,ndperiod);
        pC(i).gLGNstd=zeros(ntotal,ndperiod);
        pC(i).gEstd=zeros(ntotal,ndperiod);
        pC(i).gIstd=zeros(ntotal,ndperiod);
        pC(i).gEnstd=zeros(ntotal,ndperiod);
        pC(i).gInstd=zeros(ntotal,ndperiod);
        pC(i).gPstd=zeros(ntotal,ndperiod);
        pC(i).Vstd=zeros(ntotal,ndperiod);
        pC(i).spikes=zeros(ntotal,ndperiod);
        pC(i).gtot=zeros(ntotal,ndperiod);
        pC(i).Itot=zeros(ntotal,ndperiod);
        pC(i).Veff=zeros(ntotal,ndperiod);
        pC(i).gLGN=zeros(ntotal,ndperiod);
        pC(i).gE=zeros(ntotal,ndperiod);
        pC(i).gI=zeros(ntotal,ndperiod);
        pC(i).gEn=zeros(ntotal,ndperiod);
        pC(i).gIn=zeros(ntotal,ndperiod);
        pC(i).V=zeros(ntotal,ndperiod);
        pC(i).gP=zeros(ntotal,ndperiod);
       if i==1
        oC = pC; ipC = pC; ioC = pC;
       end
       for j=1:ntotal
       thetaPick = nP(i).indpo(j);

       pC(i).Vstd(j,:)    =  aC.Vstd(j,:,thetaPick);
       pC(i).gtotstd(j,:) =  aC.gtotstd(j,:,thetaPick);
       pC(i).Itotstd(j,:) =  aC.Itotstd(j,:,thetaPick);
       pC(i).Veffstd(j,:) =  aC.Veffstd(j,:,thetaPick);
       pC(i).gLGNstd(j,:) =  aC.gLGNstd(j,:,thetaPick);
       pC(i).gEstd(j,:)   =  aC.gEstd(j,:,thetaPick);
       pC(i).gIstd(j,:)   =  aC.gIstd(j,:,thetaPick);
       pC(i).gEnstd(j,:)  =  aC.gEnstd(j,:,thetaPick);
       pC(i).gInstd(j,:)  =  aC.gInstd(j,:,thetaPick);
       pC(i).gPstd(j,:)  =  aC.gPstd(j,:,thetaPick);
       pC(i).Vstd(j,:)  =  aC.Vstd(j,:,thetaPick);
       pC(i).spikes(j,:)  =  aC.spikes(j,:,thetaPick);
       pC(i).gtot(j,:)    =  aC.gtot(j,:,thetaPick);
       pC(i).Itot(j,:)    =  aC.Itot(j,:,thetaPick);
       pC(i).Veff(j,:)    =  aC.Veff(j,:,thetaPick);
       pC(i).gLGN(j,:)    =  aC.gLGN(j,:,thetaPick);
       pC(i).gE(j,:)      =  aC.gE(j,:,thetaPick);
       pC(i).gI(j,:)      =  aC.gI(j,:,thetaPick);
       pC(i).gEn(j,:)     =  aC.gEn(j,:,thetaPick);
       pC(i).gIn(j,:)     =  aC.gIn(j,:,thetaPick);
       pC(i).gP(j,:)     =  aC.gP(j,:,thetaPick);
       pC(i).V(j,:)       =  aC.V(j,:,thetaPick);
       pC(i).cLGN(j,:)    =  aC.cLGN(j,:,thetaPick);
       pC(i).cE(j,:)      =  aC.cE(j,:,thetaPick);
       pC(i).cI(j,:)      =  aC.cI(j,:,thetaPick);
       pC(i).cEn(j,:)     =  aC.cEn(j,:,thetaPick);
       pC(i).cIn(j,:)     =  aC.cIn(j,:,thetaPick);
       pC(i).cP(j,:)     =  aC.cP(j,:,thetaPick);
       pC(i).cLGNstd(j,:) =  aC.cLGNstd(j,:,thetaPick);
       pC(i).cEstd(j,:)   =  aC.cEstd(j,:,thetaPick);
       pC(i).cIstd(j,:)   =  aC.cIstd(j,:,thetaPick);
       pC(i).cEnstd(j,:)  =  aC.cEnstd(j,:,thetaPick);
       pC(i).cInstd(j,:)  =  aC.cInstd(j,:,thetaPick);
       pC(i).cPstd(j,:)  =  aC.cPstd(j,:,thetaPick);

       thetaPick = nP(i).indoo(j);

       oC(i).Vstd(j,:)    =  aC.Vstd(j,:,thetaPick);
       oC(i).gtotstd(j,:) =  aC.gtotstd(j,:,thetaPick);
       oC(i).Itotstd(j,:) =  aC.Itotstd(j,:,thetaPick);
       oC(i).Veffstd(j,:) =  aC.Veffstd(j,:,thetaPick);
       oC(i).gLGNstd(j,:) =  aC.gLGNstd(j,:,thetaPick);
       oC(i).gEstd(j,:)   =  aC.gEstd(j,:,thetaPick);
       oC(i).gIstd(j,:)   =  aC.gIstd(j,:,thetaPick);
       oC(i).gEnstd(j,:)  =  aC.gEnstd(j,:,thetaPick);
       oC(i).gInstd(j,:)  =  aC.gInstd(j,:,thetaPick);
       oC(i).gPstd(j,:)  =  aC.gPstd(j,:,thetaPick);
       oC(i).Vstd(j,:)  =  aC.Vstd(j,:,thetaPick);
       oC(i).spikes(j,:)  =  aC.spikes(j,:,thetaPick);
       oC(i).gtot(j,:)    =  aC.gtot(j,:,thetaPick);
       oC(i).Itot(j,:)    =  aC.Itot(j,:,thetaPick);
       oC(i).Veff(j,:)    =  aC.Veff(j,:,thetaPick);
       oC(i).gLGN(j,:)    =  aC.gLGN(j,:,thetaPick);
       oC(i).gE(j,:)      =  aC.gE(j,:,thetaPick);
       oC(i).gI(j,:)      =  aC.gI(j,:,thetaPick);
       oC(i).gEn(j,:)     =  aC.gEn(j,:,thetaPick);
       oC(i).gIn(j,:)     =  aC.gIn(j,:,thetaPick);
       oC(i).gP(j,:)     =  aC.gP(j,:,thetaPick);
       oC(i).V(j,:)       =  aC.V(j,:,thetaPick);
       oC(i).cLGN(j,:)    =  aC.cLGN(j,:,thetaPick);
       oC(i).cE(j,:)      =  aC.cE(j,:,thetaPick);
       oC(i).cI(j,:)      =  aC.cI(j,:,thetaPick);
       oC(i).cEn(j,:)     =  aC.cEn(j,:,thetaPick);
       oC(i).cIn(j,:)     =  aC.cIn(j,:,thetaPick);
       oC(i).cP(j,:)     =  aC.cP(j,:,thetaPick);
       oC(i).cLGNstd(j,:) =  aC.cLGNstd(j,:,thetaPick);
       oC(i).cEstd(j,:)   =  aC.cEstd(j,:,thetaPick);
       oC(i).cIstd(j,:)   =  aC.cIstd(j,:,thetaPick);
       oC(i).cEnstd(j,:)  =  aC.cEnstd(j,:,thetaPick);
       oC(i).cInstd(j,:)  =  aC.cInstd(j,:,thetaPick);
       oC(i).cPstd(j,:)  =  aC.cPstd(j,:,thetaPick);

       thetaPick = nP(i).indpi(j);

       ipC(i).Vstd(j,:)    =  aC.Vstd(j,:,thetaPick);
       ipC(i).gtotstd(j,:) =  aC.gtotstd(j,:,thetaPick);
       ipC(i).Itotstd(j,:) =  aC.Itotstd(j,:,thetaPick);
       ipC(i).Veffstd(j,:) =  aC.Veffstd(j,:,thetaPick);
       ipC(i).gLGNstd(j,:) =  aC.gLGNstd(j,:,thetaPick);
       ipC(i).gEstd(j,:)   =  aC.gEstd(j,:,thetaPick);
       ipC(i).gIstd(j,:)   =  aC.gIstd(j,:,thetaPick);
       ipC(i).gEnstd(j,:)  =  aC.gEnstd(j,:,thetaPick);
       ipC(i).gInstd(j,:)  =  aC.gInstd(j,:,thetaPick);
       ipC(i).gPstd(j,:)  =  aC.gPstd(j,:,thetaPick);
       ipC(i).Vstd(j,:)  =  aC.Vstd(j,:,thetaPick);
       ipC(i).spikes(j,:)  =  aC.spikes(j,:,thetaPick);
       ipC(i).gtot(j,:)    =  aC.gtot(j,:,thetaPick);
       ipC(i).Itot(j,:)    =  aC.Itot(j,:,thetaPick);
       ipC(i).Veff(j,:)    =  aC.Veff(j,:,thetaPick);
       ipC(i).gLGN(j,:)    =  aC.gLGN(j,:,thetaPick);
       ipC(i).gE(j,:)      =  aC.gE(j,:,thetaPick);
       ipC(i).gI(j,:)      =  aC.gI(j,:,thetaPick);
       ipC(i).gEn(j,:)     =  aC.gEn(j,:,thetaPick);
       ipC(i).gIn(j,:)     =  aC.gIn(j,:,thetaPick);
       ipC(i).gP(j,:)       =  aC.gP(j,:,thetaPick);
       ipC(i).V(j,:)       =  aC.V(j,:,thetaPick);
       ipC(i).cLGN(j,:)    =  aC.cLGN(j,:,thetaPick);
       ipC(i).cE(j,:)      =  aC.cE(j,:,thetaPick);
       ipC(i).cI(j,:)      =  aC.cI(j,:,thetaPick);
       ipC(i).cEn(j,:)     =  aC.cEn(j,:,thetaPick);
       ipC(i).cIn(j,:)     =  aC.cIn(j,:,thetaPick);
       ipC(i).cP(j,:)     =  aC.cP(j,:,thetaPick);
       ipC(i).cLGNstd(j,:) =  aC.cLGNstd(j,:,thetaPick);
       ipC(i).cEstd(j,:)   =  aC.cEstd(j,:,thetaPick);
       ipC(i).cIstd(j,:)   =  aC.cIstd(j,:,thetaPick);
       ipC(i).cEnstd(j,:)  =  aC.cEnstd(j,:,thetaPick);
       ipC(i).cInstd(j,:)  =  aC.cInstd(j,:,thetaPick);
       ipC(i).cPstd(j,:)  =  aC.cPstd(j,:,thetaPick);

       thetaPick = nP(i).indoi(j);

       ioC(i).Vstd(j,:)    =  aC.Vstd(j,:,thetaPick);
       ioC(i).gtotstd(j,:) =  aC.gtotstd(j,:,thetaPick);
       ioC(i).Itotstd(j,:) =  aC.Itotstd(j,:,thetaPick);
       ioC(i).Veffstd(j,:) =  aC.Veffstd(j,:,thetaPick);
       ioC(i).gLGNstd(j,:) =  aC.gLGNstd(j,:,thetaPick);
       ioC(i).gEstd(j,:)   =  aC.gEstd(j,:,thetaPick);
       ioC(i).gIstd(j,:)   =  aC.gIstd(j,:,thetaPick);
       ioC(i).gEnstd(j,:)  =  aC.gEnstd(j,:,thetaPick);
       ioC(i).gInstd(j,:)  =  aC.gInstd(j,:,thetaPick);
       ioC(i).gPstd(j,:)  =  aC.gPstd(j,:,thetaPick);
       ioC(i).Vstd(j,:)  =  aC.Vstd(j,:,thetaPick);
       ioC(i).spikes(j,:)  =  aC.spikes(j,:,thetaPick);
       ioC(i).gtot(j,:)    =  aC.gtot(j,:,thetaPick);
       ioC(i).Itot(j,:)    =  aC.Itot(j,:,thetaPick);
       ioC(i).Veff(j,:)    =  aC.Veff(j,:,thetaPick);
       ioC(i).gLGN(j,:)    =  aC.gLGN(j,:,thetaPick);
       ioC(i).gE(j,:)      =  aC.gE(j,:,thetaPick);
       ioC(i).gI(j,:)      =  aC.gI(j,:,thetaPick);
       ioC(i).gEn(j,:)     =  aC.gEn(j,:,thetaPick);
       ioC(i).gIn(j,:)     =  aC.gIn(j,:,thetaPick);
       ioC(i).gP(j,:)       =  aC.gP(j,:,thetaPick);
       ioC(i).V(j,:)       =  aC.V(j,:,thetaPick);
       ioC(i).cLGN(j,:)    =  aC.cLGN(j,:,thetaPick);
       ioC(i).cE(j,:)      =  aC.cE(j,:,thetaPick);
       ioC(i).cI(j,:)      =  aC.cI(j,:,thetaPick);
       ioC(i).cEn(j,:)     =  aC.cEn(j,:,thetaPick);
       ioC(i).cIn(j,:)     =  aC.cIn(j,:,thetaPick);
       ioC(i).cP(j,:)     =  aC.cP(j,:,thetaPick);
       ioC(i).cLGNstd(j,:) =  aC.cLGNstd(j,:,thetaPick);
       ioC(i).cEstd(j,:)   =  aC.cEstd(j,:,thetaPick);
       ioC(i).cIstd(j,:)   =  aC.cIstd(j,:,thetaPick);
       ioC(i).cEnstd(j,:)  =  aC.cEnstd(j,:,thetaPick);
       ioC(i).cInstd(j,:)  =  aC.cInstd(j,:,thetaPick);
       ioC(i).cPstd(j,:)  =  aC.cPstd(j,:,thetaPick);
        end
       clear aC;
            if (i==1)
                nP = repmat(nP,contrastLevel,1);
                tC = repmat(tC,contrastLevel,1);
                oC = repmat(oC,contrastLevel,1);
                pC = repmat(pC,contrastLevel,1);
                ioC = repmat(ioC,contrastLevel,1);
                ipC = repmat(ipC,contrastLevel,1);
            end
        [~,CVsortedID(:,i)] = sort(nP(i).cv);
        [~,PKsortedID(:,i)] = sort(nP(i).pkrate);
    end
    [~,dCVsortedID] = sort(nP(contrastLevel).cv-nP(2).cv);
    clear aC
    lgnmax_e = max(nLGN(nP(1).ei>0.5));
    lgnmin_e = min(nLGN(nP(1).ei>0.5));
    lgnmax_i = max(nLGN(nP(1).ei<0.5));
    lgnmin_i = min(nLGN(nP(1).ei<0.5));
    lowestLGN = min(lgnmin_e,lgnmin_i);
    highestLGN = max(lgnmax_e,lgnmax_i);
    [pre,preE,preI] = connections(theme,lgn,ntotal,nP(1).ei');
    toc;
    save([theme,'-tcData-x',num2str(contrastLevel),'.mat'],'tC','nP','oC','pC','ioC','ipC','conLabel','CVsortedID','PKsortedID','dCVsortedID','lgnmax_e','lgnmin_e','lgnmax_i','lgnmin_i','lowestLGN','highestLGN','pre','preE','preI');
else
    load([theme,'-tcData-x',num2str(contrastLevel),'.mat']);
end
    nv1 = size(nP(1).ei,1);
    LineScheme = {':','-.','--','-',':'};
    dtheta = 180/ntheta;
if ~fit 
    if exist([theme,'/',theme,'-fitted.mat'],'file');
        load([theme,'/',theme,'-fitted.mat']);
        disp(['fitted TC loaded']);
        tcReady = true;
    else
        tcReady = false;
    end
else
    tmpFr = zeros(p.nv1,2*ntheta,contrastLevel);
    tmpPriA = zeros(p.nv1,contrastLevel);
    for i=1:contrastLevel
        tmpFr(:,:,i) = tC(i).rate(:,1:2*ntheta);
        %tmpCV(:,i) = nP(i).cv;
        tmpPriA(:,i) = nP(i).priA;
    end
    [rmax,smax,sigfsq2] = fitTuningCurve_Quick(tmpFr,tmpPriA,contrastLevel,p.nv1,theme,ntheta,npool,theme);
    tcReady = true;
end
if tcReady
    width = sigfsq2*sqrt(log(2));
end
disp('data loaded');
    if pOSI
        meanOSI = zeros(highestLGN-lowestLGN+1,contrastLevel,2)-1;
    %    stdOSI = meanOSI;
        aeOSI = zeros(contrastLevel,nv1);
        aiOSI = zeros(contrastLevel,nv1);
        aeOSInB = aeOSI;
        aiOSInB = aiOSI;
        aiprefTheta = zeros(contrastLevel,nv1);
        aiorthTheta = aiprefTheta;
        for i = 1:contrastLevel
            aiprefTheta(i,:) = (round(nP(i).prA*180/pi/dtheta)+1)';
            aiorthTheta(i,:) = mod(aiprefTheta(i,:)+ntheta/2-1,ntheta)+1;
            %bbb = aiprefTheta(i,:)>=ntheta/2+1;
            %aiorthTheta(i,bbb) = aiprefTheta(i,bbb)-ntheta/2;
            %aiorthTheta(i,~bbb) = aiprefTheta(i,~bbb)+ntheta/2;
        end
        for i = 1:contrastLevel
            for j = 1:nv1
                nominator = tC(i).frate(j,aiprefTheta(i,j))-tC(i).frate(j,aiorthTheta(i,j));
                denorm = tC(i).frate(j,aiprefTheta(i,j))+tC(i).frate(j,aiorthTheta(i,j));
                if nP(1).ei(j) > 0.5
                    aeOSI(i,j) = nominator/denorm;
                    aeOSInB(i,j) = nominator/(denorm -2*nP(i).br(j));
                else
                    aiOSI(i,j) = nominator/denorm;
                    aiOSInB(i,j) = nominator/(denorm -2*nP(i).br(j));
                end
            end
        end
        aeOSI(isnan(aeOSI)) = 0;
        aiOSI(isnan(aiOSI)) = 0;
        aeOSInB(isnan(aeOSInB)) = 0;
        aiOSInB(isnan(aiOSInB)) = 0;

        OSIsortedID = zeros(ntotal,contrastLevel);
        for i=1:contrastLevel
            [~, OSIsortedID(:,i)] = sort([aeOSI(i,:)+aiOSI(i,:)]);
        end
        [~, dOSIsortedID] = sort((aeOSI(4,:)+aiOSI(4,:))-(aeOSI(2,:)+aiOSI(2,:)));
    end
if ~statsOnly || neuronlistOnly
    %%
    % i =32; j = 32;
    % k = (i-1)*64+j;
    % [~,k] = max(nP(1).porate);
    % k = CVsortedID(1,2);
    if isempty(neuronlist)
    %% choose over nLGN
%         level = contrastLevel;
%         fired = nP(level).pkrate(CVsortedID(:,level)) > nP(level).br(CVsortedID(:,level));
%         n = 1000;
%         sn = 3;
%         slgn = lgnmin_e:05:lgnmax_e;
%         neuronlist = zeros(sn*length(slgn)+sn,1);
%         for k = 2:length(slgn) 
%             candidates = find(nLGN(CVsortedID(:,level))>slgn(k-1) & nLGN(CVsortedID(:,level))<=slgn(k)...
%                                 & nP(level).ei(CVsortedID(:,level)) > 0.5 & fired,n,'first');
%             if ~isempty(candidates)
%                 neuronlist(sn*(k-1)+1:sn*k) = CVsortedID(candidates([1,floor(length(candidates)/2),length(candidates)]),level);
%             end
%         end
%         
%         candidates = find(nP(level).ei(CVsortedID(:,level)) < 0.5 & fired,n,'first');
%         if ~isempty(candidates)
%             neuronlist(end-sn+1:end) = CVsortedID(candidates([length(candidates),floor(length(candidates)/2),1]),level);
%         end
    %%  choose over types
        %[~,widthSortedID] = sort(sigfsq2);
        %[~,dWidthsortedID] = sort(sigfsq2(4,:)-sigfsq2(2,:));
        %sortedID = widthSortedID';
        %dsortedID = dWidthsortedID;
        sortedID = CVsortedID;
        dsortedID = dCVsortedID;
        level = contrastLevel;
        efired = nP(level).pkrate(sortedID(:,level)) > nP(level).br(sortedID(:,level)) & nP(level).ei(sortedID(:,level)) > 0.5 &...
                 nP(level).pkrate(sortedID(:,level)) > thres;
        sn = 5;
        neuronlist = zeros(sn*(p.ntypeE + p.ntypeI + 8), 1);
        ranking = cell(length(neuronlist),1);
        for i = 1:p.ntypeE
            j = (i-1)*sn;
            type0 = [p.typeE == i; false(p.nv1i,1)];
            types = type0(sortedID(:,level));
            candidates = find(efired(sortedID(:,level)) & types);
            ncand = length(candidates);
            if ncand>=sn-2
                %neuronlist(j + (1:sn-2)) = sortedID(candidates([1, floor(ncand/2), ncand]),level);
                neuronlist(j + (1:sn-2)) = sortedID(candidates([1, 2, 3]),level);
                ranking(j + (1:sn-2)) = {'smallest CV','median CV','largest CV'};
            end
            types = type0(PKsortedID(:,level));
            candidates = find(types);
            ncand = length(candidates);
            if ncand>=2
                neuronlist(j + (sn-1:sn)) = PKsortedID(candidates([floor(ncand/2),ncand]),level);
                ranking(j + (sn-1:sn)) = {'median PK','largest PK'};
            end
        end
        
        ifired = nP(level).pkrate(sortedID(:,level)) > nP(level).br(sortedID(:,level)) & nP(level).ei(sortedID(:,level)) < 0.5 &...
                 nP(level).pkrate(sortedID(:,level)) > thres;
        for i = 1:p.ntypeI
            j = p.ntypeE*sn + (i-1)*sn;
            type0 = [false(p.nv1e,1); p.typeI == i];
            types = type0(sortedID(:,4));
            candidates = find(ifired(sortedID(:,4)) & types);
            ncand = length(candidates);
            if ncand>=sn-2
                neuronlist(j + (1:sn-2)) = sortedID(candidates([1, floor(ncand/2), ncand]),level);
                ranking(j + (1:sn-2)) = {'smallest CV','median CV','largest CV'};
            end
            types = type0(PKsortedID(:,4));
            candidates = find(types);
            ncand = length(candidates);
            if ncand>=2
                neuronlist(j + (sn-1:sn)) = PKsortedID(candidates([floor(ncand/2),ncand]),level);
                ranking(j + (sn-1:sn)) = {'median PK','largest PK'};
            end
        end
            
        otherType = ' Simple';
        j = sn*(p.ntypeE+p.ntypeI);
        type0 = nP(level).sc > 1.0 & nP(level).ei > 0.5;
        types = type0(sortedID(:,level));
        candidates = find(efired(sortedID(:,level)) & types);
        ncand = length(candidates);
        if ncand>=sn-2
            %neuronlist(j + (1:sn-2)) = sortedID(candidates([1, floor(ncand/2), ncand]),level);
            neuronlist(j + (1:sn-2)) = sortedID(candidates([1, 2, 3]),level);
            ranking(j + (1:sn-2)) = strcat({'smallest CV','median CV','largest CV'},otherType);
        else
            disp(['no active exc',otherType,' is found']);
        end
        types = type0(PKsortedID(:,level));
        candidates = find(types);    
        ncand = length(candidates);
        if ncand>=2
            neuronlist(j + (sn-1:sn)) = PKsortedID(candidates([floor(ncand/2),ncand]),level);
            ranking(j + (sn-1:sn)) = strcat({'median PK','largest PK'},otherType);
        end

        otherType = ' Complex';
        j = sn*(p.ntypeE+p.ntypeI+1);
        type0 = nP(level).sc < 1.0 & nP(level).ei > 0.5; 
        types = type0(sortedID(:,level));
        candidates = find(efired(sortedID(:,level)) & types);
        ncand = length(candidates);
        if ncand>=sn-2
            %neuronlist(j + (1:sn-2)) = sortedID(candidates([1, floor(ncand/2), ncand]),level);
            neuronlist(j + (1:sn-2)) = sortedID(candidates([1, 2, 3]),level);
            ranking(j + (1:sn-2)) = strcat({'smallest CV','median CV','largest CV'},otherType);
        else
            disp(['no active exc',otherType,' is found']);
        end
        types = type0(PKsortedID(:,level));
        candidates = find(types);    
        ncand = length(candidates);
        if ncand>=2
            neuronlist(j + (sn-1:sn)) = PKsortedID(candidates([floor(ncand/2),ncand]),level);
            ranking(j + (sn-1:sn)) = strcat({'median PK','largest PK'},otherType);
        end

        otherType = ' Cort.E';
        j = sn*(p.ntypeE+p.ntypeI+2);
        ipA = round(nP(level).priA*180/pi/dtheta)+1;
        iipA = ((1:p.nv1)'-1)*(2*ntheta+1)+ipA;
        igLGN = tC(level).gLGN'.*(1+tC(level).gLGNsc');
        igE = tC(level).gE'+tC(level).gEstd';
        [~,gEsortedID] = sort(igE(iipA));
        type0 = (igLGN(iipA) < igE(iipA)) & nP(level).ei > 0.5; 
        types = type0(sortedID(:,level));
        candidates = find(efired(sortedID(:,level)) & types);
        ncand = length(candidates);
        disp([num2str(ncand),' neurons are cort.E dominated']);
        if ncand>=sn-2
            neuronlist(j + (1:sn-2)) = sortedID(candidates([1, floor(ncand/2), ncand]),level);
            ranking(j + (1:sn-2)) = strcat({'smallest CV','median CV','largest CV'},otherType);
        end
        if ncand > 0
            types = type0(gEsortedID);
            candidates = find(efired(gEsortedID) & types);    
            ncand = length(candidates);
            if ncand>=2
                neuronlist(j + (sn-1:sn)) = gEsortedID(candidates([floor(ncand/2),ncand]));
                ranking(j + (sn-1:sn)) = strcat({'median gE','largest gE'},otherType);
            else
                if ncand>0
                    neuronlist(j + (sn-1:sn)) = gEsortedID(candidates([ncand]));
                    ranking(j + (sn-1:sn)) = strcat({'the'},otherType);
                end
            end
        end

        for i = 1:p.ntypeI
            j = sn*(p.ntypeE+p.ntypeI+3);
            type0 = [false(p.nv1e,1); p.typeI == i];
            types = type0(dsortedID);
            candidates = find(ifired(dsortedID) & types); 
            ncand = length(candidates); 
            if ncand>=sn-2 
                neuronlist(j + (1:sn-2)) = dsortedID(candidates([1, floor(ncand/2), ncand]));
                ranking(j + (1:sn-2)) = {'smallest dCV','median dCV','largest dCV'};
            end
        end
        for i = 1:p.ntypeI
            j = sn*(p.ntypeE+p.ntypeI+4);
            type0 = [false(p.nv1e,1); p.typeI == i];
            types = type0(dOSIsortedID);
            candidates = find(ifired(dOSIsortedID) & types);
            ncand = length(candidates);
            if ncand>=sn-2
                neuronlist(j + (1:sn-2)) = dOSIsortedID(candidates([1, floor(ncand/2), ncand]));
                ranking(j + (1:sn-2)) = {'smallest dOSI','median dOSI','largest dOSI'};
            end
        end
        for ii = 1:2
        for i = 1:p.ntypeI
            j = sn*(p.ntypeE+p.ntypeI+4+ii);
            type0 = [false(p.nv1e,1); p.typeI == i];
            types = type0(sortedID(:,ii)); 
            candidates = find(ifired(sortedID(:,ii)) & types);
            ncand = length(candidates);
            if ncand>=sn-2
                neuronlist(j + (1:sn-2)) = sortedID(candidates([1, floor(ncand/2), ncand]),ii);
                ranking(j + (1:sn-2)) = strcat({'smallest CV','median CV','largest CV'},num2str(ii));
            end
            types = type0(PKsortedID(:,ii));
            candidates = find(types);
            ncand = length(candidates);
            if ncand>=2
                neuronlist(j + (sn-1:sn)) = PKsortedID(candidates([floor(ncand/2),ncand]),ii);
                ranking(j + (sn-1:sn)) = strcat({'median PK','largest PK'},num2str(ii));
            end
        end
        end
        for ii = 1:2
        for i = 1:p.ntypeE
            j = sn*(p.ntypeE+p.ntypeI+6+ii);
            type0 = [p.typeE==i;false(p.nv1i,1)];
            types = type0(sortedID(:,ii)); 
            candidates = find(efired(sortedID(:,ii)) & types);
            ncand = length(candidates);
            if ncand>=sn-2
                neuronlist(j + (1:sn-2)) = sortedID(candidates([1, floor(ncand/2), ncand]),ii);
                ranking(j + (1:sn-2)) = strcat({'smallest CV','median CV','largest CV'},num2str(ii));
            end
            types = type0(PKsortedID(:,ii));
            candidates = find(types);
            ncand = length(candidates);
            if ncand>=2
                neuronlist(j + (sn-1:sn)) = PKsortedID(candidates([floor(ncand/2),ncand]),ii);
                ranking(j + (sn-1:sn)) = strcat({'median PK','largest PK'},num2str(ii));
            end
        end
        end
    %%
    %%
%         candidates = find(nP(level).sc(CVsortedID(:,level))<0.6 & nP(level).ei(CVsortedID(:,level)) > 0.5,n,'first');
%         neuronlist = [CVsortedID(candidates(1:4),level); neuronlist];
    end
    pickedNeuron = neuronlist ~= 0;
    neuronlist =  neuronlist(pickedNeuron);
    ranking = ranking(pickedNeuron);
if sum(neuronlist) == 0
    disp('No single neuron match the filter');
else
    if neuronlistOnly
        return
    end
tic;
    FWHM = 2*sqrt(2*log(2));
    HWHM = sqrt(2*log(2));
    rfx = @(x0,a,b,theta,t,ra,rb) x0 + (ra.*a.*cos(t).*cos(theta)-rb.*b.*sin(t).*sin(theta));
    rfy = @(y0,a,b,theta,t,ra,rb) y0 + (ra.*a.*cos(t).*sin(theta)+rb.*b.*sin(t).*cos(theta));
    ttt = 0:0.1:2*pi;
    %finalRF(lgn,neuronlist,outputfdr);
for nn = 1:length(neuronlist)
    k = neuronlist(nn);
    
    inputAngle = nP(contrastLevel).priA(k);
    %inputAngle = nP(contrastLevel).prA(k);
    ipA = round(inputAngle*180/pi/dtheta)+1;
    if ipA < ntheta/2+1
        half = [(ipA+ntheta/2):ntheta,1:(ipA+ntheta/2)];
    else 
        half = [(ipA-ntheta/2):ntheta,1:(ipA-ntheta/2)];
    end
    if inputAngle >= pi/2
        pA = inputAngle - pi/2;
    else
        pA = inputAngle + pi/2;
    end
    grating_ipA = round(pA*180/pi/dtheta)+1;
    thetas = ((grating_ipA-ntheta/2-1):(grating_ipA+ntheta/2-1))*dtheta;

    hIndividual = figure;
    subplot(2,3,1); 
    for i=1:contrastLevel

        rho = [tC(i).frate(k,:),tC(i).frate(k,1)];
        if polarplot
            rtheta = linspace(0,2*pi,2*ntheta+1);
            polar(rtheta,rho,LineScheme{i});
        else
            plot(thetas,rho(half),'Color','k','LineStyle',LineScheme{i});
        end
        hold on
    end
    if polarplot
        polar(rtheta,nP(1).br(k)*ones(1,17),':r');
    else
        plot(thetas(end),nP(1).br(k),'*r');
        xlim([thetas(1),thetas(end)]);
    end
    ylabel('firing rate');

    if tcReady
        for il = 1:contrastLevel
            s = thetas(1):1:thetas(end);
            s0 = smax(il,k);
            if s0 < 90
                s0 = s0 + 90;
            else 
                s0 = s0 - 90;
            end
            r = rmax(il,k)*exp(-((s-s0)./sigfsq2(il,k)).^2);
            plot(s,r,LineScheme{il},'Color','g','LineWidth',2);
        end
    end
    yy = ylim();
    ylim([0,yy(2)*1.2]);
    legend(conLabel,'FontSize',FontSize-LegendOffset);

    cvL = zeros(contrastLevel,1);
    cvLnB = zeros(contrastLevel,1);
    for i=1:contrastLevel
        cvL(i) = nP(i).cv(k);
        cvLnB(i) = nP(i).cvNoBack(k);
    end
    cvLevel = strjoin(cellstr(num2str(cvL,'%1.1f'))'); 
    cvNBLevel = strjoin(cellstr(num2str(cvLnB,'%1.1f'))');
    title({['CV =', cvLevel],['evoked CV =',cvNBLevel]});

    subplot(2,3,2); hold on;
    for i=1:contrastLevel
        plot(thetas, tC(i).gLGN(k,half),'LineStyle',LineScheme{i},'Color','g');
        plot(thetas, tC(i).gLGN(k,half).*(1+tC(i).gLGNsc(k,half)),'LineStyle',LineScheme{i},'Color','k');
        plot(thetas, tC(i).gE(k,half),'LineStyle',LineScheme{i},'Color','r');
        plot(thetas, tC(i).gE(k,half)+ tC(i).gEstd(k,half),'LineStyle',LineScheme{i},'Color','b');
        plot(thetas, tC(i).gEn(k,half),'.r');
    end
    xlim([thetas(1),thetas(end)]);
    if nP(1).ei(k) > 0.5
        neuron_alias = [ranking{nn},' ',p.Etypes{p.typeE(k)}];
    else
        neuron_alias = [ranking{nn},' ',p.Itypes{p.typeI(k-p.nv1e)}];
    end
    if tcReady
        frWidth = width(:,k); 
        strFrWidth = strjoin(cellstr(num2str(frWidth,'%f'))'); 
        title({['No.', num2str(k),', ',neuron_alias],['fr HWHM = ', strFrWidth]});
    else
        title({'No.', num2str(k),', ',neuron_alias});
    end
    yy = ylim();
    ylim([0,yy(2)*1.3]);
    ylabel('Conductance');
    xlabel(['relative angle to ', conLabel{contrastLevel},' input']);
    legend({'gLGN','gLGN-F1','gE','gEstd','gEn'},'FontSize',FontSize-LegendOffset);
    
    subplot(2,3,3);hold on;
    for i=1:contrastLevel
        mmm = max(tC(i).gE(k,half));
        plot(thetas, tC(i).gE(k,half)./mmm,'LineStyle',LineScheme{i},'Color','r');
        
        mmm = max(tC(i).gI(k,half));
        plot(thetas, tC(i).gI(k,half)./mmm,'LineStyle',LineScheme{i},'Color','b');
        
        mmm = max(tC(i).gLGNsc(k,half));
        plot(thetas, tC(i).gLGNsc(k,half)./mmm,'LineStyle',LineScheme{i},'Color','k');

        mmm = max(tC(i).gLGN(k,half));
        plot(thetas, tC(i).gLGN(k,half)./mmm,'LineStyle',LineScheme{i},'Color','g');
    end
    title({['\lambda ~',num2str((nLGN(k)-lgnmin_e)/(lgnmax_e-lgnmin_e),'%1.2f'),', nLGN=',num2str(nLGN(k))],...
            ['preE: ',num2str(preE(k)),', preI: ',num2str(preI(k))]});
%     [',distance to pinwheel = ', num2str(nP(1).dist(k))];
    ylabel('Normalized g');
    xlim([thetas(1),thetas(end)]);

    subplot(2,3,4); hold on;
    for i=1:contrastLevel
        plot(thetas, tC(i).Veff(k,half),'LineStyle',LineScheme{i},'Color','b');
        plot(thetas, tC(i).Veff(k,half).*(1+tC(i).Veffsc(k,half)),'LineStyle',LineScheme{i},'Color','r');
    end
    xlim([thetas(1),thetas(end)]);
    ylabel('V_{eff}');
    yy = ylim();
    ylim([0,yy(2)*1.2]);
    legend({'Veff','Veff-F1'},'FontSize',FontSize-LegendOffset);

    OSI = zeros(contrastLevel,1);
    OSInB = OSI;
    prefinTheta = OSI;
    prefoutTheta = OSI;
    for i = 1:contrastLevel
        prefoutTheta(i) = nP(i).prA(k)*180/pi;
        iprefTheta = round(prefoutTheta(i)/dtheta+1);
        iorthTheta = mod(iprefTheta+ntheta/2-1, ntheta)+1;
        OSI(i) = (tC(i).frate(k,iprefTheta) - tC(i).frate(k,iorthTheta))/...
            (tC(i).frate(k,iprefTheta) + tC(i).frate(k,iorthTheta));
        OSInB(i) = (tC(i).frate(k,iprefTheta) - tC(i).frate(k,iorthTheta))/...
            (tC(i).frate(k,iprefTheta) + tC(i).frate(k,iorthTheta) - 2*nP(1).br(k));
        %oangles = [oangles,num2str(prefoutTheta,'%3.0f'),'^{o} '];
    end
    pick = prefoutTheta>=90;
    prefoutTheta(pick) = prefoutTheta(pick)-90;
    prefoutTheta(~pick) = prefoutTheta(~pick)+90;
    oangles = strjoin(cellstr(strcat(num2str(prefoutTheta,'%3.0f'),'^{o}'))');
    for i = 1:contrastLevel
        prefinTheta(i) = nP(i).priA(k)*180/pi;
    end
    pick = prefinTheta>=90;
    prefinTheta(pick) = prefinTheta(pick)-90;
    prefinTheta(~pick) = prefinTheta(~pick)+90;

    iangles = strjoin(cellstr(strcat(num2str(prefinTheta,'%3.0f'),'^{o}'))');

    OSI(isnan(OSI)) = 0;
    OSInB(isnan(OSInB)) = 0;
    sOSI =  strjoin(cellstr(num2str(OSI,'%1.1f'))');
    sOSInB =  strjoin(cellstr(num2str(OSInB,'%1.1f'))');
    title({['input angle :', iangles],['output angle :', oangles]});

    subplot(2,3,5); hold on;
	subregion = length(nSubLGN{k});
    subpick = 1:subregion;
    if k <= p.nv1e
        bound = p.cbounde;
        bounda = bound*ones(1,subregion);
        theta = etheta(k);
		sigma = p.esigma(k,subpick);
    	sigmb = sigma.*p.eAspectRatio(k,subpick);
        boundb = bounda.*p.eAspectRatio(k,subpick);
    else
        bound = p.cboundi;
        bounda = bound*ones(1,subregion);
        theta = itheta(k-p.nv1e);
		sigma = p.isigma(k-p.nv1e,subpick);
    	sigmb = sigma.*p.iAspectRatio(k-p.nv1e,subpick);
        boundb = bounda;
    end
    if theta >=pi/2
        gtheta = theta - pi/2;
    else
        gtheta = theta + pi/2;
    end

    peak = zeros(subregion,2);
    for j = 1:subregion
        if nSubLGN{k}(j) > 0
            peak(j,:) = mean(LGNpos(v1Map{k,j},:),1)*180/pi;
        else
            peak(j,:) = mean(LGNpos(-v1Map{k,j},:),1)*180/pi;
        end
    end
    ymax = -inf;
    ymin = inf;
    xmax = ymax;
    xmin = ymin;
    
    for j = 1:subregion
        xxx0 = rfx(peak(j,1)/180*pi,sigmb(j),sigma(j),gtheta,ttt,boundb(j),bounda(j))*180/pi;         
        yyy0 = rfy(peak(j,2)/180*pi,sigmb(j),sigma(j),gtheta,ttt,boundb(j),bounda(j))*180/pi;
        xxx1 = rfx(peak(j,1)/180*pi,sigmb(j),sigma(j),gtheta,ttt,HWHM,HWHM)*180/pi;         
        yyy1 = rfy(peak(j,2)/180*pi,sigmb(j),sigma(j),gtheta,ttt,HWHM,HWHM)*180/pi;
        xxx2 = rfx(peak(j,1)/180*pi,sigmb(j),sigma(j),gtheta,ttt,FWHM,FWHM)*180/pi;         
        yyy2 = rfy(peak(j,2)/180*pi,sigmb(j),sigma(j),gtheta,ttt,FWHM,FWHM)*180/pi;
        ymin = min(min(yyy2),ymin);
        ymax = max(max(yyy2),ymax);
        xmin = min(min(xxx2),xmin);
        xmax = max(max(xxx2),xmax);
        x = [xmin,xmax];
        if nSubLGN{k}(j) > 0   
            plot(xxx2,yyy2,'--r','LineWidth',1);
            plot(xxx1,yyy1,':r','LineWidth',1); 
            plot(xxx0,yyy0,'-r','LineWidth',1); 
            plot(peak(j,1),peak(j,2),'sr');
            %plot3(LGNpos(v1Map{k,j},1)*180/pi, LGNpos(v1Map{k,j},2)*180/pi,lgnStrength{k,j}*5,'*r','MarkerSize',4.0);
            plot(LGNpos(v1Map{k,j},1)*180/pi, LGNpos(v1Map{k,j},2)*180/pi,'*r','MarkerSize',4.0);
            plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':r');
        else
            plot(xxx2,yyy2,'--b','LineWidth',1);
            plot(xxx1,yyy1,':b','LineWidth',1); 
            plot(xxx0,yyy0,'-b','LineWidth',1); 
            plot(peak(j,1),peak(j,2),'sb');
            %plot3(LGNpos(-v1Map{k,j},1)*180/pi, LGNpos(-v1Map{k,j},2)*180/pi,lgnStrength{k,j}*5,'ob','MarkerSize',4.0);
            plot(LGNpos(-v1Map{k,j},1)*180/pi, LGNpos(-v1Map{k,j},2)*180/pi,'ob','MarkerSize',4.0);
            plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':b');
        end
    end
    grid on
    axis equal
    ylim([ymin,ymax]);
    xlim([xmin,xmax]);
    dxy = p.rlgn*180/pi;
    edge = 1;
    x = (xmin-edge*dxy):(dxy/20):(xmax+edge*dxy);
    y = (ymin-edge*dxy):(dxy/20):(ymax+edge*dxy);   

%     zmax = 0;
    for j = 1:subregion
       if nSubLGN{k}(j) > 0
            [X, Y, Z] = multiLGNspatialKernel(LGNpos(v1Map{k,j},1)*180/pi, LGNpos(v1Map{k,j},2)*180/pi,x,y,nSubLGN{k}(j),lgnStrength{k,j},200);
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','m','LineWidth',2);
            newZmax = max(max(Z));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','m','LineWidth',2);
       else
            [X, Y, Z] = multiLGNspatialKernel(LGNpos(-v1Map{k,j},1)*180/pi, LGNpos(-v1Map{k,j},2)*180/pi,x,y,-nSubLGN{k}(j),lgnStrength{k,j},200);
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','g','LineWidth',2);
            newZmax = max(max(Z));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','g','LineWidth',2);
       end
%        mesh(X,Y,Z,'FaceColor','none');
%        if  newZmax > zmax
%            zmax = newZmax;
%        end
    end
    ylim([y(1),y(end)]);
    xlim([x(1),x(end)]);
    xlabel('degree');
    ylabel('degree');
    %campos([xrange(1),yrange(1),zmax*1.5]);

    scL = zeros(contrastLevel,1);
    for i = 1:contrastLevel
        scL(i) = nP(i).sc(k);
    end
    scLevel = strjoin(cellstr(num2str(scL,'%1.1f'))');

    title({['preset grating angle: ',num2str(gtheta*180/pi,'%3.0f'),'^{o}, #sub:',num2str(reshape(nSubLGN{k},[1,subregion]),'%2i')],...
			['F1/F0 =', scLevel]});

    subplot(2,3,6); hold on;
    for i=1:contrastLevel
        errorbar(thetas+5*(i-1), tC(i).gI(k,half),tC(i).gIstd(k,half),'LineStyle',LineScheme{i},'Color','b');
        plot(thetas, tC(i).gIn(k,half),'.b');
    end
    title({['OSI: ',sOSI],['OSInB: ',sOSInB]});
    xlim([thetas(1),thetas(end)]);
    yy = ylim();
    ylim([0,yy(2)*1.2]);
    legend({'gI','gIn'},'FontSize',FontSize-LegendOffset);
    if ~isempty(format)
        set(gcf,'Renderer','Painters')
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
        if strcmp(format,'fig')
            saveas(hIndividual,[outputfdr,'/',num2str(k),'-',neuron_alias,'-tc-',theme,'.',format]);
        else
            print(hIndividual,[outputfdr,'/',num2str(k),'-',neuron_alias,'-tc-',theme,'.',format],printDriver,dpi);
        end
    end
    %% phase
    %if k>p.nv1e
    oiLabel = {'o','i'};
    offset = 0;
    if Phase
        for oi = 1:2
            if ~operiod && oi == 1
                continue
            end
            if oi == 1
                iitheta = nP(level).indpo(k);
            else
                iitheta = nP(level).indpi(k);
            end
            level = contrastLevel;
            hPhase = figure;
            %colormap(hPhase,phaseSkip(ndperiod+1));
            %subplot(2,2,1)
            %x = 1;
            %y = p.ev1y;
            c = zeros(p.nv1e,1);
            frC = zeros(preE(k),1)-1;
            for i = 1:preE(k)
                id = pre(k).ID_exc(i);
                titheta = nP(level).indpo(id);
                if abs(titheta - iitheta) <= offset || abs(titheta - iitheta) >= ntheta - offset
                    [frC(i),c(id)] = max(pC(level).spikes(id,:));
                end
            end
            %c=reshape(c,[p.ev1y,p.ev1x]);
            c = c./ndperiod;
            ch = c(c>0).*ndperiod;
            frC = frC.*ndperiod/rt.*(profiles(pre(k).s_exc))';
            frC = frC(frC>=0);
            %imagesc(x,y,c);
            %title(['Exc-',num2str((iitheta-1)*dtheta),'^o']);
           
            if ~isempty(ch) 
             
            hExc = subplot(2,2,1);
            nfr = 10;
            ctrs = cell(2,1);
            ctrs{1} = 1:ndperiod;
            lctrsx = length(ctrs{1});
            lctrsy = nfr+1;
            dTickX = 5/ndperiod;
            dTickY = 0.2;
            tickPosX = 0.5:lctrsx*dTickX:lctrsx+0.5;
            tickPosY = 0.5:lctrsy*dTickY:lctrsy+0.5;
            tickLabelX = num2str((linspace(0,2*pi,length(tickPosX)))'.*180/pi);

            maxfr = max(frC);
            minfr = min(frC);
            ctrs{2} = linspace(minfr,maxfr,lctrsy);
            tickLabelY = flipud(num2str((linspace(minfr,maxfr,length(tickPosY)))'));
            pair = [ch,frC];
            den = hist3(pair,ctrs);
            den = den(1:ntheta,1:nfr);
            den = den/max(max(den));
            imagesc([1,lctrsx],[lctrsy,1],den');
            set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
            xlabel('Phase');
            ylabel('Excitation');
            title([num2str(k),'-',neuron_alias,' Dist. of phase and Excitation',num2str(iitheta)]);
            colormap(hExc,redOnly);
            else
            title(['no iso-oriented exc neuron were connected']);
            end
            
            c = zeros(p.nv1i,1);
            frC = zeros(preI(k),1)-1;
            for i = 1:preI(k)
                id = pre(k).ID_inh(i);
                titheta = nP(level).indpo(id);
                if abs(titheta - iitheta) <= offset || abs(titheta - iitheta) >= ntheta - offset
                    [frC(i),c(id)] = max(pC(level).spikes(id,:));
                end
            end
            %c=reshape(c,[p.iv1y,p.iv1x]);
            c = c./ndperiod;
            ch = c(c>0).*ndperiod;
            frC = frC.*ndperiod/rt;
            frC = frC(frC>=0);

            if ~isempty(ch) 
            hInh=subplot(2,2,2);
            maxfr = max(frC);
            minfr = min(frC);
            ctrs{2} = linspace(minfr,maxfr,lctrsy);
            tickLabelY = flipud(num2str((linspace(minfr,maxfr,length(tickPosY)))'));
            pair = [ch,frC];
            den = hist3(pair,ctrs);
            den = den(1:ntheta,1:nfr);
            den = den/max(max(den));
            imagesc([1,lctrsx],[lctrsy,1],den');
            set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
            xlabel('Phase');
            ylabel('Hz');
            title([num2str(k),'-',neuron_alias,' Dist. of phase and Inh FR',num2str(iitheta)]);
            colormap(hInh,blueOnly);
            else
            title(['no iso-oriented inh neuron were connected']);
            end

            %%
            dddtheta = pi/ntheta;
            nfr = 20;
            ctrs = cell(2,1);
            ctrs{1} = 0:dddtheta:pi;
            lctrsx = length(ctrs{1});
            lctrsy = nfr+1;
            dTickX = 2/ntheta;
            dTickY = 0.2;
            tickPosX = 0.5:lctrsx*dTickX:lctrsx+0.5;
            tickPosY = 0.5:lctrsy*dTickY:lctrsy+0.5;
            tickLabelX = num2str((linspace(0,pi,length(tickPosX)))'.*180/pi);

            hExc = subplot(2,2,3);
            frTarget = tC(level).frate(pre(k).ID_exc,iitheta);
            thetaTarget = nP(level).prA(pre(k).ID_exc);
            maxfr = max(frTarget);
            minfr = min(frTarget);
            ctrs{2} = linspace(minfr,maxfr,lctrsy);
            tickLabelY = flipud(num2str((linspace(minfr,maxfr,length(tickPosY)))'));
            pair = [thetaTarget,frTarget];
            den = hist3(pair,ctrs);
            den = den(1:ntheta,1:nfr);
            den = den/max(max(den));
            imagesc([1,lctrsx],[lctrsy,1],den');
            set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
            xlabel('\theta');
            ylabel('FR Hz');
            colormap(hExc,redOnly);

            hInh = subplot(2,2,4);
            frTarget = tC(level).frate(pre(k).ID_inh,iitheta);
            thetaTarget = nP(level).prA(pre(k).ID_inh);
            maxfr = max(frTarget);
            minfr = min(frTarget);
            ctrs{2} = linspace(minfr,maxfr,lctrsy);
            tickLabelY = flipud(num2str((linspace(minfr,maxfr,length(tickPosY)))'));
            pair = [thetaTarget,frTarget];
            den = hist3(pair,ctrs);
            den = den(1:ntheta,1:nfr);
            den = den/max(max(den));
            imagesc([1,lctrsx],[lctrsy,1],den');
            set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
            xlabel('\theta');
            ylabel('FR Hz');
            colormap(hInh,blueOnly);

            if ~isempty(format)
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                if strcmp(format,'fig')
                    saveas(hPhase,[outputfdr,'/',num2str(k),'-',neuron_alias,'-C',num2str(level),'-FRphaseMap',oiLabel{oi},'-A',num2str(iitheta),'-',theme,'.',format]);
                else
                    print(hPhase,[outputfdr,'/',num2str(k),'-',neuron_alias,'-C',num2str(level),'-FRphaseMap',oiLabel{oi},'-A',num2str(iitheta),'-',theme,'.',format],printDriver,dpi);
                end
            end
        end
    end
    %end
    %% current
    if current
        hCurrent = figure;
        for i=1:contrastLevel
            subplot(2,2,1)
            hold on
            errorbar(thetas+5*(i-1), tC(i).currEt(k,half),tC(i).currEtstd(k,half),'LineStyle',LineScheme{i},'Color','g');
            if i == 1, title('LGN current');xlim([thetas(1),thetas(end)]);end
            
            subplot(2,2,2)
            hold on
            errorbar(thetas+5*(i-1), tC(i).currEn(k,half),tC(i).currEnstd(k,half),'LineStyle',LineScheme{i},'Color','r');
            errorbar(thetas+5*(i-1), tC(i).currIn(k,half),tC(i).currInstd(k,half),'LineStyle',LineScheme{i},'Color','b');
            if i == 1, title('E and I noise current');xlim([thetas(1),thetas(end)]);end
            
            subplot(2,2,3)
            hold on
            errorbar(thetas+5*(i-1), tC(i).currEc(k,half),tC(i).currEcstd(k,half),'LineStyle',LineScheme{i},'Color','r');
            if i == 1, title('Exc Cortical current');xlim([thetas(1),thetas(end)]);end
            
            subplot(2,2,4)
            hold on
            errorbar(thetas+5*(i-1), tC(i).currIc(k,half),tC(i).currIcstd(k,half),'LineStyle',LineScheme{i},'Color','b');
            if i == 1, title('Inh Cortical current');xlim([thetas(1),thetas(end)]);end
        end
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(hCurrent,[outputfdr,'/',num2str(k),'-',neuron_alias,'-current-',theme,'.',format]);
            else
                print(hCurrent,[outputfdr,'/',num2str(k),'-',neuron_alias,'-current-',theme,'.',format],printDriver,dpi);
            end
        end
    end
    %% timeline average
    if meanTimeLine
        if operiod
            hTimeline = figure;
            for i=1:contrastLevel
                
                iprefTheta = round((nP(i).prA(k)*180/pi)/dtheta+1);
                subplot(contrastLevel,4,(i-1)*4+1);hold on;
                
                errorbar(1:ndperiod,pC(i).gLGN(k,:),pC(i).gLGNstd(k,:),'-g');
                plot(ndperiod+0.6,mean(pC(i).gLGN(k,:)),'sg');
                errorbar((1:ndperiod)+0.3,oC(i).gLGN(k,:),oC(i).gLGNstd(k,:),':k');
                plot(ndperiod+0.3,mean(oC(i).gLGN(k,:)),'^k');
                if i == 1
                    title(['gLGN, F1/F0: ',num2str(tC(i).gLGNsc(k,iprefTheta),'%1.2f')]); 
                else
                    title(num2str(tC(i).gLGNsc(k,iprefTheta),'%1.2f')); 
                end
                ylabel(conLabel{i});
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+2);hold on;
                
                errorbar(1:ndperiod,pC(i).gI(k,:), pC(i).gIstd(k,:),'-b');
                plot(ndperiod+0.6,mean(pC(i).gI(k,:)),'sb');
                errorbar(1:ndperiod,oC(i).gI(k,:), oC(i).gIstd(k,:),':k');
                plot(ndperiod,mean(oC(i).gI(k,:)),'^k');
                if i == 1, title(['gInh, ',neuron_alias]); end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+3); hold on;  
                errorbar(1:ndperiod,pC(i).gE(k,:), pC(i).gEstd(k,:),'-r');
                plot(ndperiod+0.6,mean(pC(i).gE(k,:)),'sr');
                errorbar((1:ndperiod)+0.3,oC(i).gE(k,:), oC(i).gEstd(k,:),':k');
                plot(ndperiod+0.3,mean(oC(i).gE(k,:)),'^k');
                if i == 1
                    title(['gE, F1/F0: ',num2str(tC(i).gEsc(k,iprefTheta),'%1.2f')]);
                else
                    title(num2str(tC(i).gEsc(k,iprefTheta),'%1.2f'));
                end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+4); hold on;
                if i==1,title(['FR, ID.', num2str(k)]);end
                plot(1:ndperiod,pC(i).spikes(k,:)*ndperiod/rt,'-k');
                plot(1:ndperiod,oC(i).spikes(k,:)*ndperiod/rt,':k');
                plot(ndperiod+0.6,mean(pC(i).spikes(k,:))*ndperiod/rt,'sk');
                plot(ndperiod+0.3,mean(oC(i).spikes(k,:))*ndperiod/rt,'^k');
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                if i==contrastLevel, xlabel('phase');end
            end
            if ~isempty(format)
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                if strcmp(format,'fig')
                    saveas(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-operiod-',theme,'.',format]);
                else
                    print(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-operiod-',theme,'.',format],printDriver,dpi);
                end
            end
        end

        hTimeline = figure;
        for i=1:contrastLevel
            
            iprefTheta = round((nP(i).priA(k)*180/pi)/dtheta+1);
            subplot(contrastLevel,4,(i-1)*4+1);hold on;
            
            errorbar(1:ndperiod,ipC(i).gLGN(k,:),ipC(i).gLGNstd(k,:),'-g');
            plot(ndperiod+0.6,mean(ipC(i).gLGN(k,:)),'sg');
            errorbar((1:ndperiod)+0.3,ioC(i).gLGN(k,:),ioC(i).gLGNstd(k,:),':k');
            plot(ndperiod+0.3,mean(ioC(i).gLGN(k,:)),'^k');
            if i == 1
                title(['gLGN, F1/F0: ',num2str(tC(i).gLGNsc(k,iprefTheta),'%1.2f')]); 
            else
                title(num2str(tC(i).gLGNsc(k,iprefTheta),'%1.2f')); 
            end
            ylabel(conLabel{i});
            ylim([0,inf]);
            xlim([0,ndperiod+2]);
            
            subplot(contrastLevel,4,(i-1)*4+2);hold on;
            
            errorbar(1:ndperiod,ipC(i).gI(k,:), ipC(i).gIstd(k,:),'-b');
            plot(ndperiod+0.6,mean(ipC(i).gI(k,:)),'sb');
            errorbar(1:ndperiod,ioC(i).gI(k,:), ioC(i).gIstd(k,:),':k');
            plot(ndperiod,mean(ioC(i).gI(k,:)),'^k');
            if i == 1, title(['gInh, ',neuron_alias]); end
            ylim([0,inf]);
            xlim([0,ndperiod+2]);
            
            subplot(contrastLevel,4,(i-1)*4+3); hold on;  
            errorbar(1:ndperiod,ipC(i).gE(k,:), ipC(i).gEstd(k,:),'-r');
            plot(ndperiod+0.6,mean(ipC(i).gE(k,:)),'sr');
            errorbar((1:ndperiod)+0.3,ioC(i).gE(k,:), ioC(i).gEstd(k,:),':k');
            plot(ndperiod+0.3,mean(ioC(i).gE(k,:)),'^k');
            if i == 1
                title(['gE, F1/F0: ',num2str(tC(i).gEsc(k,iprefTheta),'%1.2f')]);
            else
                title(num2str(tC(i).gEsc(k,iprefTheta),'%1.2f'));
            end
            ylim([0,inf]);
            xlim([0,ndperiod+2]);
            
            subplot(contrastLevel,4,(i-1)*4+4); hold on;
            if i==1,title(['FR, ID.', num2str(k)]);end
            plot(1:ndperiod,ipC(i).spikes(k,:)*ndperiod/rt,'-k');
            plot(1:ndperiod,ioC(i).spikes(k,:)*ndperiod/rt,':k');
            plot(ndperiod+0.6,mean(ipC(i).spikes(k,:))*ndperiod/rt,'sk');
            plot(ndperiod+0.3,mean(ioC(i).spikes(k,:))*ndperiod/rt,'^k');
            ylim([0,inf]);
            xlim([0,ndperiod+2]);
            if i==contrastLevel, xlabel('phase');end
        end
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-iperiod-',theme,'.',format]);
            else
                print(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-iperiod-',theme,'.',format],printDriver,dpi);
            end
        end
        if current
            if operiod
                hTimeline = figure;
                for i=1:contrastLevel
                    
                    iprefTheta = round((nP(i).prA(k)*180/pi)/dtheta+1);
                    subplot(contrastLevel,4,(i-1)*4+1);hold on;
                    
                    errorbar(1:ndperiod,pC(i).cLGN(k,:),pC(i).cLGNstd(k,:),'-g');
                    plot(ndperiod+0.6,mean(pC(i).cLGN(k,:)),'sg');
                    errorbar((1:ndperiod)+0.3,oC(i).cLGN(k,:),oC(i).cLGNstd(k,:),':k');
                    plot(ndperiod+0.3,mean(oC(i).cLGN(k,:)),'^k');
                    if i == 1
                        title('cLGN'); 
                    end
                    ylabel(conLabel{i});
                    ylim([0,inf]);
                    xlim([0,ndperiod+2]);
                    
                    subplot(contrastLevel,4,(i-1)*4+2);hold on;
                    
                    errorbar(1:ndperiod,pC(i).cIstd(k,:), pC(i).cIstd(k,:),'-b');
                    plot(ndperiod+0.6,mean(pC(i).cIstd(k,:)),'sb');
                    errorbar(1:ndperiod,oC(i).cI(k,:), oC(i).cIstd(k,:),':k');
                    plot(ndperiod,mean(oC(i).cI(k,:)),'^k');
                    if i == 1, title(['cInh, ',neuron_alias]); end
                    ylim([0,inf]);
                    xlim([0,ndperiod+2]);
                    
                    subplot(contrastLevel,4,(i-1)*4+3); hold on;  
                    errorbar(1:ndperiod,pC(i).cE(k,:), pC(i).cEstd(k,:),'-r');
                    plot(ndperiod+0.6,mean(pC(i).cE(k,:)),'sr');
                    errorbar((1:ndperiod)+0.3,oC(i).cE(k,:), oC(i).cEstd(k,:),':k');
                    plot(ndperiod+0.3,mean(oC(i).cE(k,:)),'^k');
                    if i == 1
                        title('cExc');
                    end
                    ylim([0,inf]);
                    xlim([0,ndperiod+2]);
                    
                    subplot(contrastLevel,4,(i-1)*4+4); hold on;
                    if i==1,title(['FR, ID.', num2str(k)]);end
                    plot(1:ndperiod,pC(i).spikes(k,:)*ndperiod/rt,'-k');
                    plot(1:ndperiod,oC(i).spikes(k,:)*ndperiod/rt,':k');
                    plot(ndperiod+0.6,mean(pC(i).spikes(k,:))*ndperiod/rt,'sk');
                    plot(ndperiod+0.3,mean(oC(i).spikes(k,:))*ndperiod/rt,'^k');
                    ylim([0,inf]);
                    xlim([0,ndperiod+2]);
                    if i==contrastLevel, xlabel('phase');end
                end
                if ~isempty(format)
                    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                    if strcmp(format,'fig')
                        saveas(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-coperiod-',theme,'.',format]);
                    else
                        print(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-coperiod-',theme,'.',format],printDriver,dpi);
                    end
                end
            end
            hTimeline = figure;
            for i=1:contrastLevel
                
                iprefTheta = round((nP(i).prA(k)*180/pi)/dtheta+1);
                subplot(contrastLevel,4,(i-1)*4+1);hold on;
                
                errorbar(1:ndperiod,ipC(i).cLGN(k,:),ipC(i).cLGNstd(k,:),'-g');
                plot(ndperiod+0.6,mean(ipC(i).cLGN(k,:)),'sg');
                errorbar((1:ndperiod)+0.3,ioC(i).cLGN(k,:),ioC(i).cLGNstd(k,:),':k');
                plot(ndperiod+0.3,mean(ioC(i).cLGN(k,:)),'^k');
                if i == 1
                    title('cLGN'); 
                end
                ylabel(conLabel{i});
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+2);hold on;
                
                errorbar(1:ndperiod,ipC(i).cIstd(k,:), ipC(i).cIstd(k,:),'-b');
                plot(ndperiod+0.6,mean(ipC(i).cIstd(k,:)),'sb');
                errorbar(1:ndperiod,ioC(i).cI(k,:), ioC(i).cIstd(k,:),':k');
                plot(ndperiod,mean(ioC(i).cI(k,:)),'^k');
                if i == 1, title(['cInh, ',neuron_alias]); end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+3); hold on;  
                errorbar(1:ndperiod,ipC(i).cE(k,:), ipC(i).cEstd(k,:),'-r');
                plot(ndperiod+0.6,mean(ipC(i).cE(k,:)),'sr');
                errorbar((1:ndperiod)+0.3,ioC(i).cE(k,:), ioC(i).cEstd(k,:),':k');
                plot(ndperiod+0.3,mean(ioC(i).cE(k,:)),'^k');
                if i == 1
                    title('cExc');
                end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+4); hold on;
                if i==1,title(['FR, ID.', num2str(k)]);end
                plot(1:ndperiod,ipC(i).spikes(k,:)*ndperiod/rt,'-k');
                plot(1:ndperiod,ioC(i).spikes(k,:)*ndperiod/rt,':k');
                plot(ndperiod+0.6,mean(ipC(i).spikes(k,:))*ndperiod/rt,'sk');
                plot(ndperiod+0.3,mean(ioC(i).spikes(k,:))*ndperiod/rt,'^k');
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                if i==contrastLevel, xlabel('phase');end
            end
            if ~isempty(format)
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                if strcmp(format,'fig')
                    saveas(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-ciperiod-',theme,'.',format]);
                else
                    print(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-ciperiod-',theme,'.',format],printDriver,dpi);
                end
            end
        end
    end
    if dThetaSize
        for oi = 1:2
            if ~operiod && oi == 1
                continue;
            end
            hdThetaSizeFR = figure; 
            ndtheta = ntheta/2;
            dThetaColorScheme = ones(ndtheta,3);
            dThetaColorScheme(:,1) = (linspace(0,11/15,ndtheta))';
            dThetaColorScheme(:,2) = ones(ndtheta,1);
            dThetaColorScheme(:,3) = ones(ndtheta,1);
            for i=1:contrastLevel
                if oi == 1
                    iiitheta = nP(i).indpo(k);
                else
                    iiitheta = nP(i).indpi(k);
                end
                %subplot(contrastLevel,2,2*(i-1)+1)
                subplot(contrastLevel,1,i)
                hold on
                dThetaVec = abs(nP(i).prA(pre(k).ID_exc)-nP(i).priA(k));

                tempPick = dThetaVec > pi/2;
                dThetaVec(tempPick) = pi - dThetaVec(tempPick);
                idThetaVec = (round(dThetaVec*180/pi/dtheta)+1)';
                FRvec = tC(i).frate(pre(k).ID_exc,iiitheta);
                SizeVec = pre(k).s_exc;
                for j=1:ndtheta
                    target = idThetaVec==j-1;
                    plot(SizeVec(target),FRvec(target),'o','Color',hsv2rgb(dThetaColorScheme(j,:)));
                end
                xlabel('PSPSize'); ylabel('FR');title('dTheta 0->90, red->blue');
            end

            if ~isempty(format)
                    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                    if strcmp(format,'fig')
                        saveas(hdThetaSizeFR,[outputfdr,'/',num2str(k),'-',neuron_alias,'-C',num2str(level),'-dThetaSizeFR',oiLabel{oi},'-',num2str(iitheta),'-',theme,'.',format]);
                    else
                        print(hdThetaSizeFR,[outputfdr,'/',num2str(k),'-',neuron_alias,'-C',num2str(level),'-dThetaSizeFR',oiLabel{oi},'-',num2str(iitheta),'-',theme,'.',format],printDriver,dpi);
                    end
            end
        end
    end
end
disp('individual examples done');
toc;
end
end
%%
if (~individOnly)
tic;
%hPhase = figure;
%colormap(hPhase,phaseNSkip(ndperiod+1));
%level = contrastLevel;
%subplot(1,2,1)
%x = 1;
%y = p.ev1y;
%c = zeros(p.nv1e,1)-1;
%for i=1:p.ev1x
%    for j=1:p.ev1y
%        id = j+(i-1)*p.ev1y;
%        [~,c(id)] = max(pC(level).spikes(id,:));
%    end
%end
%c=reshape(c,[p.ev1y,p.ev1x]);
%c = c./ndperiod;
%imagesc(x,y,c);
%title('Exc');
%
%subplot(1,2,2)
%x = 1;
%y = p.iv1y;
%c = zeros(p.nv1i,1)-1;
%for i=1:p.iv1x
%    for j=1:p.iv1y
%        id = j+(i-1)*p.iv1y;
%        [~,c(id)] = max(pC(level).spikes(p.nv1e+id,:));
%    end
%end
%c = reshape(c,[p.iv1y,p.iv1x]);
%c = c./ndperiod;
%imagesc(x,y,c);
%title('Inh');
%if ~isempty(format)
%    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
%    print(hPhase,[outputfdr,'/','C',num2str(level),'-FRphaseMap-',theme,'.',format],printDriver,dpi);
%end
halfNdperiod = (ndperiod + mod(ndperiod,2))/2;
pe = 1:p.nv1e;

hPhase = figure;
tauCorr = zeros(p.nv1e,ndperiod);
edges = (0:ndperiod)*360/ndperiod;
for i =1:contrastLevel
    subplot(contrastLevel,2,(i-1)*2+1)
        meanEtemp = mean(ipC(i).cE(pe,:),2)*ones(1,ndperiod);
        meanItemp = mean(ipC(i).cI(pe,:),2)*ones(1,ndperiod);
        varEtemp = var(ipC(i).cE(pe,:),1,2);
        varItemp = var(ipC(i).cI(pe,:),1,2);
        for j = 1:ndperiod
            tauCorr(:,j) = mean((ipC(i).cE(pe,:)-meanEtemp).*circshift(ipC(i).cI(pe,:)-meanItemp,j-1,2),2)./sqrt(varEtemp.*varItemp);
        end
        [~,idphase] = max(tauCorr,[],2);
        %pt = idphase > halfNdperiod;
        %idphase(pt) = ndperiod - idphase(pt);
        histogram(idphase*360/ndperiod,edges);
        if i==contrastLevel, xlabel('\DeltaPhase');end
        if i==1, title('Phase maxCorr(E_{input},I_{input})');end
        ylabel('# Neurons');
    subplot(contrastLevel,2,(i-1)*2+2);
        meanEtemp = mean(ipC(i).cLGN(pe,:),2)*ones(1,ndperiod);
        varEtemp = var(ipC(i).cLGN(pe,:),1,2);
        for j = 1:ndperiod
            tauCorr(:,j) = mean((ipC(i).cLGN(pe,:)-meanEtemp).*circshift(ipC(i).cI(pe,:)-meanItemp,j-1,2),2)./sqrt(varEtemp.*varItemp);
        end
        [~,idphase] = max(tauCorr,[],2);
        %pt = idphase > halfNdperiod;
        %idphase(pt) = ndperiod - idphase(pt);
        histogram(idphase*360/ndperiod,edges);
        if i==contrastLevel, xlabel('\DeltaPhase');end
        if i==1, title('(LGN_{input},I_{input})');end
end

if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hPhase,[outputfdr,'/popInputPhase-',theme,'.',format]);
    else
        print(hPhase,[outputfdr,'/popInputPhase-',theme,'.',format],printDriver,dpi);
    end
end

hInput = figure;
alignedPhase = (halfNdperiod + mod(halfNdperiod,2))/2;
for i =1:contrastLevel
    tLGN = zeros(ndperiod,2); 
    tE = zeros(ndperiod,2); 
    tI = zeros(ndperiod,2);
    tP = zeros(ndperiod,2); 
    tEn = zeros(ndperiod,2); 
    tIn = zeros(ndperiod,2); 
    tL = zeros(ndperiod,1);
    subplot(contrastLevel,3,(i-1)*3+1)
    hold on
        [~,imax] = max(ipC(i).cLGN(pe,:),[],2);
        for j = 1:p.nv1e
            shift = alignedPhase-imax(j);
            tLGN(:,1) = tLGN(:,1) - circshift(ipC(i).cLGN(j,:),shift,2)';
            tLGN(:,2) = tLGN(:,2) + circshift(ipC(i).cLGNstd(j,:),shift,2)';
            tE(:,1) = tE(:,1) - circshift(ipC(i).cE(j,:),shift,2)';
            tE(:,2) = tE(:,2) + circshift(ipC(i).cEstd(j,:),shift,2)';
            tI(:,1) = tI(:,1) - circshift(ipC(i).cI(j,:),shift,2)';
            tI(:,2) = tI(:,2) + circshift(ipC(i).cIstd(j,:),shift,2)';
            tP(:,1) = tP(:,1) - circshift(ipC(i).cP(j,:),shift,2)';
            tP(:,2) = tP(:,2) + circshift(ipC(i).cPstd(j,:),shift,2)';
            tEn(:,1) = tEn(:,1) - circshift(ipC(i).cEn(j,:),shift,2)';
            tEn(:,2) = tEn(:,2) + circshift(ipC(i).cEnstd(j,:),shift,2)';
            tIn(:,1) = tIn(:,1) - circshift(ipC(i).cIn(j,:),shift,2)';
            tIn(:,2) = tIn(:,2) + circshift(ipC(i).cInstd(j,:),shift,2)';
            tL = tL - gL*(circshift(ipC(i).V(j,:),shift,2)-vrest)';
        end
        tL = tL./p.nv1e; tLGN = tLGN./p.nv1e; tE = tE./p.nv1e; tI = tI./p.nv1e; tEn = tEn./p.nv1e; tIn = tIn./p.nv1e; tP= tP./p.nv1e;
        x = (0:ndperiod-1).*360/ndperiod;
        plot(x,tL,'.','MarkerSize',3);
        errorbar(x,tLGN(:,1),tLGN(:,2),'-g');
        errorbar(x,tE(:,1),tE(:,2),'-r');
        errorbar(x,tI(:,1),tI(:,2),'-b');
        errorbar(x,tEn(:,1),tEn(:,2),':m','LineWidth',0.5);
        errorbar(x,tIn(:,1),tIn(:,2),':c','LineWidth',0.5);
        errorbar(x,tP(:,1),tP(:,2),':k');
        plot(x,tLGN(:,1)+tE(:,1)+tI(:,1)+tEn(:,1)+tIn(:,1)+tP(:,1)+tL,'-k');
        xlim([0,360]);
        xlabel('Phase (deg)');
        ylabel('Current');

    subplot(contrastLevel,3,(i-1)*3+2)
    hold on
        target = ipC(i).cLGN(pe,:) + ipC(i).cE(pe,:) + ipC(i).cI(pe,:) + ipC(i).cP(pe,:) + ipC(i).cEn(pe,:) + ipC(i).cIn(pe,:) + gL*(ipC(i).V(pe,:)-vrest);
        [Allcounts,Alledges] = histcounts(mean(-target),'Normalization','pdf');
        [Ecounts,Eedges] = histcounts(mean(-ipC(i).cE(pe,:)),'Normalization','pdf');
        [Icounts,Iedges] = histcounts(mean(-ipC(i).cI(pe,:)),'Normalization','pdf');
        [Pcounts,Pedges] = histcounts(mean(-ipC(i).cP(pe,:)),'Normalization','pdf');
        [LGNcounts,LGNedges] = histcounts(mean(-ipC(i).cLGN(pe,:)),'Normalization','pdf');
        plot((Eedges(1:end-1)+Eedges(2:end))./2,Ecounts,'r');
        plot((Iedges(1:end-1)+Iedges(2:end))./2,Icounts,'b');
        plot((Pedges(1:end-1)+Pedges(2:end))./2,Pcounts,':k');
        plot((Alledges(1:end-1)+Alledges(2:end))./2,Allcounts,'k');
        plot((LGNedges(1:end-1)+LGNedges(2:end))./2,LGNcounts,'g');
        xlabel('avg. Current');
        ylabel('% Neurons');
        if i==1, title('Pref. Input');end

    subplot(contrastLevel,3,(i-1)*3+3)
    hold on
        target = ioC(i).cLGN(pe,:) + ioC(i).cE(pe,:) + ioC(i).cI(pe,:) + ioC(i).cP(pe,:) + ioC(i).cEn(pe,:) + ioC(i).cIn(pe,:) + gL*(ioC(i).V(pe,:)-vrest);
        [Allcounts,Alledges] = histcounts(mean(-target),'Normalization','pdf');
        [Ecounts,Eedges] = histcounts(mean(-ioC(i).cE(pe,:)),'Normalization','pdf');
        [Icounts,Iedges] = histcounts(mean(-ioC(i).cI(pe,:)),'Normalization','pdf');
        [Pcounts,Pedges] = histcounts(mean(-ioC(i).cP(pe,:)),'Normalization','pdf');
        [LGNcounts,LGNedges] = histcounts(mean(-ioC(i).cLGN(pe,:)),'Normalization','pdf');
        plot((Eedges(1:end-1)+Eedges(2:end))./2,Ecounts,'-r');
        plot((Iedges(1:end-1)+Iedges(2:end))./2,Icounts,'-b');
        plot((Pedges(1:end-1)+Pedges(2:end))./2,Pcounts,':k');
        plot((Alledges(1:end-1)+Alledges(2:end))./2,Allcounts,'-k');
        plot((LGNedges(1:end-1)+LGNedges(2:end))./2,LGNcounts,'-g');
        xlabel('avg. Current');
        ylabel('% Neurons');
        if i==1, title('Orth. Input');end
end
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hInput,[outputfdr,'/Input-',theme,'.',format]);
    else
        print(hInput,[outputfdr,'/Input-',theme,'.',format],printDriver,dpi);
    end
end

hgExc = figure;
nd = p.enormDistance;
lnd = 16;
edges = linspace(0,1,lnd);
xxx = (edges(1:lnd-1) + edges(2:lnd))/2;
[~,~,ind] = histcounts(nd,edges);
gExc = zeros(lnd-1,2);
frExc = zeros(lnd-1,2);
for i=1:lnd-1
    target = pC(contrastLevel).gE(1:p.nv1e,:);
    target = target(ind==i,:);
    target = target';
    gExc(i,1) = mean(target(:));
    gExc(i,2) = mean(std(target));
    target = pC(contrastLevel).spikes(1:p.nv1e,:)*ndperiod/rt;
    target = target(ind==i,:);
    target = target';
    frExc(i,1) = mean(target(:));
    frExc(i,2) = mean(std(target));
end
subplot(2,1,1)
errorbar(xxx,gExc(:,1),gExc(:,2));
ylabel('gE');
subplot(2,1,2)
errorbar(xxx,frExc(:,1),frExc(:,2));
ylabel('frt');
xlabel('norm''Distance');

if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hgExc,[outputfdr,'/gExc-',theme,'.',format]);
    else
        print(hgExc,[outputfdr,'/gExc-',theme,'.',format],printDriver,dpi);
    end
end
ll = ((1:ntheta)-1)/ntheta*pi;
ll = ones(p.nv1,1)*ll;
for i=1:contrastLevel
    target = tC(i).gLGN + tC(i).gE;
    nP(i).giCV    = 1-abs(sum(tC(i).gI(:,1:ntheta).*exp(2i*ll),2)./sum(tC(i).gI(:,1:ntheta),2));
    nP(i).geCV    = 1-abs(sum(tC(i).gE(:,1:ntheta).*exp(2i*ll),2)./sum(tC(i).gE(:,1:ntheta),2));
    nP(i).gepCV   = 1-abs(sum(tC(i).gLGN(:,1:ntheta).*exp(2i*ll),2)./sum(tC(i).gLGN(:,1:ntheta),2));
    nP(i).getotCV = 1-abs(sum(target(:,1:ntheta).*exp(2i*ll),2)./sum(target(:,1:ntheta),2));
end
hgCVheat = figure;

p1 =contrastLevel-2; p2 = contrastLevel;
while p1<=0
    p1 = p1 + 1;
end
nCV = 8;
ctrs = cell(2,1);
dTick = 0.2;
pick = nP(contrastLevel).ei>0.5; 

hExc = subplot(2,2,1);
pair1 = 1-nP(p1).geCV(pick);
pair2 = 1-nP(p2).geCV(pick);
minCtrs = min(min(pair1),min(pair2));
maxCtrs = max(max(pair1),max(pair2));
if minCtrs~=maxCtrs && ~isnan(minCtrs) && ~isnan(maxCtrs)
[tick,n0] = autoAxis(minCtrs,maxCtrs,round(1/dTick),[0,1]);
tickLabel = num2str(tick');
lctrsx = (n0-1)*nCV+1;
lctrsy = lctrsx;
ctrs{1} = linspace(tick(1),tick(end),lctrsx);
ctrs{2} = ctrs{1};
tickPosY = 0.5:nCV:(lctrsy-1+0.5);
tickPosX = 0.5:nCV:(lctrsx-1+0.5);

CVpair = [pair1, pair2];
denCVpair = hist3(CVpair,ctrs);
maxDen = max(max(denCVpair));
denCVpair = denCVpair/maxDen;
denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
hold on
plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);

set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
title('gEcortical');
colormap(hExc,redOnly);
end

hEp = subplot(2,2,2);
pair1 = 1-nP(p1).gepCV(pick);
pair2 = 1-nP(p2).gepCV(pick);
minCtrs = min(min(pair1),min(pair2));
maxCtrs = max(max(pair1),max(pair2));
if minCtrs~=maxCtrs && ~isnan(minCtrs) && ~isnan(maxCtrs)
[tick,n0] = autoAxis(minCtrs,maxCtrs,round(1/dTick),[0,1]);
tickLabel = num2str(tick');
lctrsx = (n0-1)*nCV+1;
lctrsy = lctrsx;
ctrs{1} = linspace(tick(1),tick(end),lctrsx);
ctrs{2} = ctrs{1};
tickPosY = 0.5:nCV:(lctrsx-1+0.5);
tickPosX = 0.5:nCV:(lctrsy-1+0.5);
CVpair = [pair1, pair2];
denCVpair = hist3(CVpair,ctrs);
maxDen = max(max(denCVpair));
denCVpair = denCVpair/maxDen;
denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
hold on
plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);

set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
title('gLGN');
colormap(hEp,redOnly);
end

hEtot = subplot(2,2,3);
pair1 = 1-nP(p1).getotCV(pick);
pair2 = 1-nP(p2).getotCV(pick);
minCtrs = min(min(pair1),min(pair2));
maxCtrs = max(max(pair1),max(pair2));
if minCtrs~=maxCtrs && ~isnan(minCtrs) && ~isnan(maxCtrs)
[tick,n0] = autoAxis(minCtrs,maxCtrs,round(1/dTick),[0,1]);
tickLabel = num2str(tick');
lctrsx = (n0-1)*nCV+1;
lctrsy = lctrsx;
ctrs{1} = linspace(tick(1),tick(end),lctrsx);
ctrs{2} = ctrs{1};
tickPosY = 0.5:nCV:(lctrsx-1+0.5);
tickPosX = 0.5:nCV:(lctrsy-1+0.5);
CVpair = [pair1, pair2];
denCVpair = hist3(CVpair,ctrs);
maxDen = max(max(denCVpair));
denCVpair = denCVpair/maxDen;
denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
hold on
plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);

set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
title('gEtot');
colormap(hEtot,redOnly);
end

hInh = subplot(2,2,4);
pick = nP(contrastLevel).ei>0.5;
pair1 = 1-nP(p1).giCV(pick);
pair2 = 1-nP(p2).giCV(pick);
minCtrs = min(min(pair1),min(pair2));
maxCtrs = max(max(pair1),max(pair2));
if minCtrs~=maxCtrs && ~isnan(minCtrs) && ~isnan(maxCtrs)
[tick,n0] = autoAxis(minCtrs,maxCtrs,round(1/dTick),[0,1]);
tickLabel = num2str(tick');
lctrsx = (n0-1)*nCV+1;
lctrsy = lctrsx;
ctrs{1} = linspace(tick(1),tick(end),lctrsx);
ctrs{2} = ctrs{1};
tickPosY = 0.5:nCV:(lctrsx-1+0.5);
tickPosX = 0.5:nCV:(lctrsy-1+0.5);

CVpair = [pair1, pair2];
denCVpair = hist3(CVpair,ctrs);
maxDen = max(max(denCVpair));
denCVpair = denCVpair/maxDen;
denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
hold on
plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);

set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
title('gInh');
colormap(hInh,blueOnly);
end
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hgCVheat,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_gCVheat-',theme,'.',format]);
    else
        print(hgCVheat,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_gCVheat-',theme,'.',format],printDriver,dpi);
    end
end

hCVpair = figure;

cvLGNmean = zeros(highestLGN-lowestLGN+1,contrastLevel,2)-1;
%cvLGNstd = cvLGNmean;
p1 =contrastLevel-2; p2 = contrastLevel;
while p1<=0
    p1 = p1 + 1;
end
%suptitle(theme);
subplot(1,2,1); hold on;
for i=lgnmin_e:lgnmax_e
    pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    if ~isempty(pick)
        plot(1-nP(p1).cv(pick),1-nP(p2).cv(pick),'o','Color',[0.1+0.899*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
    end
    for j = 1:contrastLevel
        pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(j).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            cvLGNmean(i-lowestLGN+1,j,1) = mean(nP(j).cv(pick));
    %        cvLGNstd(i-lowestLGN+1,j,1) = std(nP(j).cv(pick));
        end
    end
end
plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
ylabel(['1-CV at ',num2str(12.5*2^(p2-1),'%3.1f'),'% Contrast']);
xlabel(['1-CV at ',num2str(12.5*2^(p1-1),'%3.1f'),'% Contrast']);
title('Excitatory');
subplot(1,2,2); hold on;

for i=lgnmin_i:lgnmax_i
    pick= nP(contrastLevel).ei<0.5 & nLGN==i & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(contrastLevel).pkrate > thres;
    if ~isempty(pick)
        plot(1-nP(p1).cv(pick),1-nP(p2).cv(pick),'o','Color',[0,0,0.1+0.899*(i-lgnmin_i)/(lgnmax_i-lgnmin_i)]);
    end
    for j = 1:contrastLevel
        pick= nP(contrastLevel).ei<0.5 & nLGN==i & nP(j).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            cvLGNmean(i-lowestLGN+1,j,2) = mean(nP(j).cv(pick));
%            cvLGNstd(i-lowestLGN+1,j,2) = std(nP(j).cv(pick));
        end
    end
end
plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
ylabel(['1-CV at ',conLabel{p2},' Contrast']);
xlabel(['1-CV at ',conLabel{p1},' Contrast']);
title('Inhibitory');
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hCVpair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_CV-',theme,'.',format]);
    else
        print(hCVpair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_CV-',theme,'.',format],printDriver,dpi);
    end
end
%%
hCVheatPair = figure;

p1 =contrastLevel-2; p2 = contrastLevel;
while p1<=0
    p1 = p1 + 1;
end
nCV = 8;
ctrs = cell(2,1);
dTick = 0.2;
hExc = subplot(1,2,1);
pick = nP(contrastLevel).ei>0.5 & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(p1).pkrate>thres & nP(p2).pkrate>thres;
if sum(pick) > 0
    pair1 = 1-nP(p1).cv(pick);
    pair2 = 1-nP(p2).cv(pick);
    minCtrs = min(min(pair1),min(pair2));
    maxCtrs = max(max(pair1),max(pair2));
    [tick,n0] = autoAxis(minCtrs,maxCtrs,round(1/dTick),[0,1]);
    tickLabel = num2str(tick');
    lctrsx = (n0-1)*nCV+1;
    lctrsy = lctrsx;
    ctrs{1} = linspace(tick(1),tick(end),lctrsx);
    ctrs{2} = ctrs{1};
    tickPosY = 0.5:nCV:(lctrsy-1+0.5);
    tickPosX = 0.5:nCV:(lctrsx-1+0.5);
    
    CVpair = [pair1, pair2];
    denCVpair = hist3(CVpair,ctrs);
    maxDen = max(max(denCVpair));
    denCVpair = denCVpair/maxDen;
    denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
    imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
    hold on
    plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
    title('Excitatory');
    xlabel('gOSI(25%)')
    ylabel('gOSI(100%)')
    daspect([1,1,1]);
    
    set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
    colormap(hExc,redOnly);
end

hInh = subplot(1,2,2);
pick = nP(contrastLevel).ei<0.5 & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(p1).pkrate>thres & nP(p2).pkrate>thres;
if sum(pick) >0
    pair1 = 1-nP(p1).cv(pick);
    pair2 = 1-nP(p2).cv(pick);
    minCtrs = min(min(pair1),min(pair2));
    maxCtrs = max(max(pair1),max(pair2));
    [tick,n0] = autoAxis(minCtrs,maxCtrs,round(1/dTick),[0,1]);
    tickLabel = num2str(tick');
    nCV = 5;
    lctrsx = (n0-1)*nCV+1;
    lctrsy = lctrsx;
    ctrs{1} = linspace(tick(1),tick(end),lctrsx);
    ctrs{2} = ctrs{1};
    tickPosY = 0.5:nCV:(lctrsx-1+0.5);
    tickPosX = 0.5:nCV:(lctrsy-1+0.5);
    
    CVpair = [pair1, pair2];
    denCVpair = hist3(CVpair,ctrs);
    maxDen = max(max(denCVpair));
    denCVpair = denCVpair/maxDen;
    denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
    imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
    hold on
    plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
    %plot(0.5:lctrsx+0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
    title('Inhibitory');
    xlabel('gOSI(25%)')
    ylabel('gOSI(100%)')
    daspect([1,1,1]);
    
    set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
    colormap(hInh,blueOnly);
end

if ~isempty(format)
    %set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hCVheatPair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_CVheat-',theme,'.',format]);
    else
        print(hCVheatPair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_CVheat-',theme,'.',format],printDriver,dpi);
    end
end
hCVnLGN = figure;
subplot(1,2,1); 
hold on;
% bar(0:lgnmax,cvLGNmean(:,:,1));
x = (-0.2:(highestLGN-lowestLGN-0.2));
for i=1:contrastLevel
%    x = (-0.2:(highestLGN-lowestLGN-0.2))+0.1*(i-1);
    pick = cvLGNmean(:,i,1)>=0;
    plot(x(pick),1-cvLGNmean(pick,i,1),'-*','Color',[0.1+0.899*i/contrastLevel,0,0]);
end
xlabel('# LGN input');
ylabel('1-CV');
title('active Exc 1-CV');
xlim([lgnmin_e, lgnmax_e]);
ylim([0,inf]);
yy = ylim();
ylim([yy(1),yy(2)*1.2]);
legend(conLabel);

subplot(1,2,2); 
hold on;
for i=1:contrastLevel
    pick = cvLGNmean(:,i,2)>=0;
    plot(x(pick),1-cvLGNmean(pick,i,2),'-*','Color',[0,0,0.1+0.899*i/contrastLevel]);
end
xlabel('# LGN input');
ylabel('1-CV');
title('active Inh 1-CV');
xlim([lgnmin_i, lgnmax_i]);
ylim([0,inf]);
yy = ylim();
ylim([yy(1),yy(2)*1.2]);
legend(conLabel);

if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hCVnLGN,[outputfdr,'/','CVoverLGN-',theme,'.',format]);
    else
        print(hCVnLGN,[outputfdr,'/','CVoverLGN-',theme,'.',format],printDriver,dpi);
    end
end
% subplot(1,2,2); hold on;
% for i=1:contrastLevel
%     x = 0:lgnmax+0.2;
%     errorbar(x,cvLGNmean(:,i,2),cvLGNstd(:,i,2),'*','Color',[0,0,0.1+0.9*i/lgnmax]);
%     bar(x,cvLGNmean(i,:,2),'FaceColor',[0.1+0.9*i/lgnmax,0,0],'BarWidth',0.1);
% end
%% OSI
hPrefOrth = figure; 
eoorate = zeros(contrastLevel,p.nv1e);
eporate = zeros(contrastLevel,p.nv1e);
eoirate = zeros(contrastLevel,p.nv1e);
epirate = zeros(contrastLevel,p.nv1e);

ioorate = zeros(contrastLevel,p.nv1i);
iporate = zeros(contrastLevel,p.nv1i);
ioirate = zeros(contrastLevel,p.nv1i);
ipirate = zeros(contrastLevel,p.nv1i);
for j=1:contrastLevel
    nntheta = size(tC(j).frate,2);
    %pref
    pick = false(ntotal * nntheta,1);
    opPick = aiprefTheta(j,:) + (0:nntheta:((ntotal-1)*nntheta));
    pick(opPick) = true;
    pick = reshape(pick,[nntheta,ntotal])';
    
    eipick = pick;
    eipick(nP(1).ei < 0.5,:) = false;
    eporate(j,:) = tC(j).frate(eipick)';
    
    eipick = pick;
    eipick(nP(1).ei > 0.5,:) = false;
    iporate(j,:) = tC(j).frate(eipick)';
    % orth
    pick = false(ntotal * nntheta,1);
    opPick = aiorthTheta(j,:) + (0:nntheta:((ntotal-1)*nntheta));
    pick(opPick) = true;
    pick = reshape(pick,[nntheta,ntotal])';
    
    eipick = pick;
    eipick(nP(1).ei < 0.5,:) = false;
    eoorate(j,:) = tC(j).frate(eipick)';
    
    eipick = pick;
    eipick(nP(1).ei > 0.5,:) = false;
    ioorate(j,:) = tC(j).frate(eipick)';
    
    clear eipick opPick
end


aiprefTheta_in = zeros(contrastLevel,nv1);
aiorthTheta_in = aiprefTheta_in;
for i = 1:contrastLevel
    aiprefTheta_in(i,:) = (round(nP(i).priA*180/pi/dtheta)+1)';
    aiorthTheta_in(i,:) = mod(aiprefTheta_in(i,:)+ntheta/2-1,ntheta)+1;
    %bbb = aiprefTheta(i,:)>=ntheta/2+1;
    %aiorthTheta(i,bbb) = aiprefTheta(i,bbb)-ntheta/2;
    %aiorthTheta(i,~bbb) = aiprefTheta(i,~bbb)+ntheta/2;
end

for j=1:contrastLevel
    nntheta = size(tC(j).frate,2);
    %pref
    pick = false(ntotal * nntheta,1);
    opPick = aiprefTheta_in(j,:) + (0:nntheta:((ntotal-1)*nntheta));
    pick(opPick) = true;
    pick = reshape(pick,[nntheta,ntotal])';
    
    eipick = pick;
    eipick(nP(1).ei < 0.5,:) = false;
    epirate(j,:) = tC(j).frate(eipick)';
    
    eipick = pick;
    eipick(nP(1).ei > 0.5,:) = false;
    ipirate(j,:) = tC(j).frate(eipick)';
    % orth
    pick = false(ntotal * nntheta,1);
    opPick = aiorthTheta_in(j,:) + (0:nntheta:((ntotal-1)*nntheta));
    pick(opPick) = true;
    pick = reshape(pick,[nntheta,ntotal])';
    
    eipick = pick;
    eipick(nP(1).ei < 0.5,:) = false;
    eoirate(j,:) = tC(j).frate(eipick)';
    
    eipick = pick;
    eipick(nP(1).ei > 0.5,:) = false;
    ioirate(j,:) = tC(j).frate(eipick)';
    
    clear eipick opPick 
end


cc = contrastLevel;
while cc <= 0
    cc = cc + 1;
end

shifted = 4;
subplot(2,2,1); hold on
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eporate,2),std(eporate,1,2),'o:r');
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoorate,2),std(eoorate,1,2),'o:k');
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(epirate,2),std(epirate,1,2),'o:m');
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoirate,2),std(eoirate,1,2),'o:g');
legend({'Prefout','Orthout','Prefin','Orthin'});

epick = nP(1).ei > 0.5;
pick = nP(cc).pkrate(epick) > nP(1).br(epick) & nP(contrastLevel).pkrate(epick) > thres;
npick = sum(pick);
if npick > 0
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'^-r');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'^-k');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'^-m');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'^-g');
end
ylim([0,inf]);xlim([0,110]);
%legend('Prefered Orientation','Orthogonal Orientation')
xlabel('Contrast Level %');
ylabel('Firing Rate');
title(['Exc ',num2str(npick/p.nv1e*100,'%3.1f'),'% evoked']);

subplot(2,2,2); hold on
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(iporate,2),std(iporate,1,2),'o:r');
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(ioorate,2),std(ioorate,1,2),'o:k');
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(ipirate,2),std(ipirate,1,2),'o:m');
errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(ioirate,2),std(ioirate,1,2),'o:g');

ipick = nP(1).ei < 0.5;
pick = nP(cc).pkrate(ipick) > nP(1).br(ipick) & nP(contrastLevel).pkrate(ipick) > thres;
npick = sum(pick);
if npick > 0
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(iporate(:,pick),2),std(iporate(:,pick),1,2),'^-r');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(ioorate(:,pick),2),std(ioorate(:,pick),1,2),'^-k');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(ipirate(:,pick),2),std(ipirate(:,pick),1,2),'^-m');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(ioirate(:,pick),2),std(ioirate(:,pick),1,2),'^-g');
end
ylim([0,inf]);
%legend('Prefered Orientation','Orthogonal Orientation')
xlabel('Contrast Level');
ylabel('Firing Rate');
title(['Inh ',num2str(npick/p.nv1i*100,'%3.1f'),'% evoked']);

subplot(2,2,3); hold on

epick = nP(cc).ei>0.5;
pick = nP(cc).sc(epick)>1.0;
npickAll = sum(pick);
if npickAll > 0
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'o:r');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'o:k');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'o:m');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'o:g');
end

pick = nP(cc).sc(epick)>1.0 & nP(cc).pkrate(epick) > nP(1).br(epick) & nP(contrastLevel).pkrate(epick) > thres;
npick = sum(pick);
if npick > 0
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'^-r');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'^-k');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'^-m');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'^-g');
end

ylim([0,inf]);xlim([0,110]);
%legend('Prefered Orientation','Orthogonal Orientation')
xlabel('Contrast Level %');
ylabel('Firing Rate');
title(['Simple Exc at ', conLabel{cc},' Contrast ', num2str(npick/npickAll*100,'%3.1f'),'% evoked']);

subplot(2,2,4); hold on
pick = nP(cc).sc(epick)<1.0;
npickAll = sum(pick);
if npickAll > 0
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'o:r');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'o:k');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'o:m');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'o:g');
end

pick = nP(cc).sc(epick)<1.0 & nP(cc).pkrate(epick) > nP(1).br(epick) & nP(contrastLevel).pkrate(epick) > thres;
npick = sum(pick);
disp([num2str(npick),' active complex cells at ',conLabel{cc},' contrast ', num2str(npickAll),' non-active complex cells']);
if npick > 0
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'^-r');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'^-k');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'^-m');
    errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'^-g');
end

ylim([0,inf]);xlim([0,110]);
%legend('Prefered Orientation','Orthogonal Orientation')
xlabel('Contrast Level %');
ylabel('Firing Rate');
title(['Complex Exc at ', conLabel{cc},' Contrast ', num2str(npick/npickAll*100,'%3.1f'),'% evoked']);

if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hPrefOrth,[outputfdr,'/','GainCurve-',theme,'.',format]);
    else
        print(hPrefOrth,[outputfdr,'/','GainCurve-',theme,'.',format],printDriver,dpi);
    end
end

clear eporate eoorate ioorate iporate
if pOSI
    hOSIpair = figure;
    subplot(1,2,1); hold on;
    p1 =contrastLevel-2; p2 = contrastLevel;
    while p1<=0
        p1 = p1+1;
    end
    %suptitle(theme);
    for i=lgnmin_e:lgnmax_e
        pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            plot(aeOSI(p1,pick),aeOSI(p2,pick),'o','Color',[0.1+0.899*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
        end
        for j = 1:contrastLevel
            pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(j).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
            if ~isempty(pick)
                meanOSI(i-lowestLGN+1,j,1) = mean(aeOSI(j,pick));
            end
        end
    end
    plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
    
    xlabel(['OSI at ',conLabel{p2},' Contrast']);
    ylabel(['OSI at ',conLabel{p1},' Contrast']);

    title('Excitatory');
    subplot(1,2,2); hold on;
    for i=lgnmin_i:lgnmax_i
        pick= nP(contrastLevel).ei<0.5 & nLGN==i & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            plot(aiOSI(p1,pick),aiOSI(p2,pick),'o','Color',[0,0,0.1+0.899*(i-lgnmin_i)/(lgnmax_i-lgnmin_i)]);
        end
        for j = 1:contrastLevel
            pick= nP(contrastLevel).ei<0.5 & nLGN==i & nP(j).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
            if ~isempty(pick)
                meanOSI(i-lowestLGN+1,j,2) = mean(aiOSI(j,pick));
            end
        end
    end
    plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
    xlabel(['OSI at ',conLabel{p2},' Contrast']);
    ylabel(['OSI at ',conLabel{p1},' Contrast']);
    title('Inhibitory');

%    subplot(2,2,3); hold on;
%    pick= nP(contrastLevel).ei>0.5 & nP(p1).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
%    binCountsE = hist(aeOSI(p1,nP(1).ei>0.5),0.05:0.1:0.95);
%    pick= nP(contrastLevel).ei<0.5 & nP(p1).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
%    binCountsI = hist(aiOSI(p1,nP(1).ei<0.5),0.05:0.1:0.95);
%
%    bar(0.05:0.1:0.95,[binCountsI;binCountsE]'./ntotal.*100);
%    legend({'Inh','Exc'});
%    ylabel('% neuron');
%    xlabel(['OSI at ',num2str(12.5*2^(p2-1),'%3.1f'),'% Contrast']);
%    subplot(2,2,4); hold on;
%    title(['meanOSI = ',num2str(mean(aeOSI(4,nP(4).sc<1.0)),'%1.1f')]);
%    hist(aeOSI(4,nP(4).sc<1.0),0.05:0.1:0.95);
%    ylabel('# Complex neurons');
%    xlabel('OSI');
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hOSIpair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_OSI-',theme,'.',format]);
        else
            print(hOSIpair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_OSI-',theme,'.',format],printDriver,dpi);
        end
    end
    %% OSI over # LGN
    hOSIovernLGN_OSIvsCV = figure; 

    subplot(2,2,1); hold on;
    x = (-0.2:(highestLGN-lowestLGN-0.2));
    for i=1:contrastLevel
        %x = (-0.2:(highestLGN-lowestLGN-0.2))+0.1*(i-1);
        pick = meanOSI(:,i,1)>=0;
        plot(x(pick),meanOSI(pick,i,1),'-*','Color',[0.1+0.899*i/contrastLevel,0,0]);
    end
    xlabel('# LGN input'); ylabel('mean OSI');
    title('active Exc OSI');
    xlim([lgnmin_e, lgnmax_e]);
    ylim([0,inf]);
    yy = ylim();
    ylim([yy(1),yy(2)*1.2]);
    legend(conLabel);

    subplot(2,2,3); hold on;
    for i=1:contrastLevel
        pick = meanOSI(:,i,2)>=0;
        plot(x(pick),meanOSI(pick,i,2),'-*','Color',[0,0,0.1+0.899*i/contrastLevel]);
    end
    xlabel('# LGN input'); ylabel('mean OSI');
    title('active Inh OSI');
    xlim([lgnmin_i, lgnmax_i]);
    ylim([0,inf]);
    yy = ylim();
    ylim([yy(1),yy(2)*1.2]);
    legend(conLabel);

    subplot(1,2,2); hold on;
    pp = 3;
    for i=lgnmin_e:lgnmax_e
        pick= nP(1).ei>0.5 & nLGN==i & nP(pp).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        plot(aeOSI(3,pick),1-nP(pp).cv(pick),'o','Color',[0.1+0.899*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
    end
    for i=lgnmin_i:lgnmax_i
        pick= nP(1).ei<0.5 & nLGN==i & nP(pp).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        plot(aiOSI(3,pick),1-nP(pp).cv(pick),'o','Color',[0,0,0.1+0.899*(i-lgnmin_i)/(lgnmax_i-lgnmin_i)]);
    end

    plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
    xlabel('OSI');
    ylabel('1-CV');
    title(['1-CV vs. OSI at',num2str(12.5*2^(pp-1),'%3.1f'),'% Contrast']);
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hOSIovernLGN_OSIvsCV,[outputfdr,'/','OSIovernLGN_OSIvsCV-',theme,'.',format]);
        else
            print(hOSIovernLGN_OSIvsCV,[outputfdr,'/','OSIovernLGN_OSIvsCV-',theme,'.',format],printDriver,dpi);
        end
    end
    clear aeOSI aiOSI aeOSInB aiOSInB
end
%% simple complex
hSimCom = figure;
meanSC = zeros(highestLGN-lowestLGN+1,contrastLevel,2)-1;

p1 =contrastLevel-2; p2 = contrastLevel;
while p1<=0
    p1 = p1+1;
end
subplot(1,2,1); hold on;
for i=lgnmin_e:lgnmax_e
    pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    if ~isempty(pick)
        plot(nP(p1).sc(pick),nP(p2).sc(pick),'o','Color',[0.1+0.899*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
    end
    for j = 1:contrastLevel
        pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(j).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            meanSC(i-lowestLGN+1,j,1) = mean(nP(j).sc(pick));
        end
    end
end
x = 0:0.1:2;
nx = length(x);
plot(x,linspace(0,2,nx),'-.k','LineWidth',2);
plot(x,ones(nx,1),'-.k','LineWidth',2);
plot(ones(nx,1),x,'-.k','LineWidth',2);

xlabel(['F1/F0 at ',conLabel{p2},' Contrast']);
ylabel(['F1/F0 at ',conLabel{p1},' Contrast']);

title('Excitatory');
subplot(1,2,2); hold on;
for i=lgnmin_i:lgnmax_i
    pick= nP(contrastLevel).ei<0.5 & nLGN==i & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    if ~isempty(pick)
        plot(nP(p1).sc(pick),nP(p2).sc(pick),'o','Color',[0,0,0.1+0.899*(i-lgnmin_i)/(lgnmax_i-lgnmin_i)]);
    end
    for j = 1:contrastLevel
        pick= nP(contrastLevel).ei<0.5 & nLGN==i & nP(j).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            meanSC(i-lowestLGN+1,j,2) = mean(nP(j).sc(pick));
        end
    end
end
plot(x,linspace(0,2,nx),'-.k','LineWidth',2);
plot(x,ones(nx,1),'-.k','LineWidth',2);
plot(ones(nx,1),x,'-.k','LineWidth',2);

xlabel(['F1/F0 at ',conLabel{p2},' Contrast']);
ylabel(['F1/F0 at ',conLabel{p1},' Contrast']);
title('Inhibitory');
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hSimCom,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_F1overF2-',theme,'.',format]);
    else
        print(hSimCom,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_F1overF2-',theme,'.',format],printDriver,dpi);
    end
end
clear meanSC
%% orientations
hOriExc = figure;

thetaRange = 180/ntheta*((0:ntheta)-0.5);
for j = 1:contrastLevel
    pick = nP(contrastLevel).ei>0.5 & nP(j).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    oriRate = zeros(ntheta,1);
    po1d = nP(j).prA(pick)*180/pi;
    subplot(contrastLevel,2,(j-1)*2+1);
    [binCount, ind] = histc(po1d,thetaRange);
    binCount = binCount(1:ntheta);
    fired_pkrate = nP(j).pkrate(pick);
    for i=1:ntheta
        thetaPick = (ind == i);
        if sum(thetaPick)>0
            oriRate(i) = mean(fired_pkrate(thetaPick));
        end
    end
    plot((0:(ntheta-1))*180/ntheta,oriRate,'LineWidth',2);
    xlim([0,180-180/ntheta]);
    xlabel('pref angle (^{o})');
    ylabel('Hz');

    subplot(contrastLevel,2,(j-1)*2+2);
    plot((0:(ntheta-1))*180/ntheta,binCount,'LineWidth',2);
    xlim([0,180-180/ntheta]);
    xlabel('pref angle (^{o})');
    ylabel('# Exc Cell');
    title([num2str(sum(pick)/p.nv1e*100,'%2.1f'),'% Exc evoked']);
end

if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hOriExc,[outputfdr,'/','ExcOrientation-',theme,'.',format]);
    else
        print(hOriExc,[outputfdr,'/','ExcOrientation-',theme,'.',format],printDriver,dpi);
    end
end
hOriCmp = figure; hold on;
p1 = 3; p2 = 4;
shifts = zeros(lgnmax_e-lgnmin_e+1,1);
subplot(1,2,1); hold on
% amax = zeros(contrastLevel,nv1);
if tcReady 
    amax = smax;
%else
%    amax = (aiprefTheta-1)*dtheta; 
%end
    amax(amax>180) = amax(amax>180)-180;
    amax(amax<0) = amax(amax<0)+180;
    for i=lgnmin_e:lgnmax_e
        pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(contrastLevel).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            plot(amax(p1,pick),amax(p2,pick),'o','Color',[0.1+0.899*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
            shifts(i-lgnmin_e+1) = sum(abs(amax(p2,pick)-amax(p1,pick)))/sum(pick);
        end
    end
    xlabel(['Preferred Orientation at ',num2str(12.5*2^(p1-1),'%3.1f'),'% Contrast (degree)']); 
    ylabel(['Preferred Orientation at ',num2str(12.5*2^(p2-1),'%3.1f'),'% Contrast']);
    xlim([0,210]);
else
    pick= nP(contrastLevel).ei>0.5 & nP(contrastLevel).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    [y,x]=hist((aiprefTheta(contrastLevel,:)-1)*dtheta,-dtheta:dtheta:180-dtheta);
    bar(x,y);
    [y,x]=hist((aiprefTheta(contrastLevel,pick)-1)*dtheta,-dtheta:dtheta:180-dtheta);
    hBar = bar(x,y);
    set(hBar,'FaceColor','r');
    legend({'Inh','Exc'});
    xlabel('Preferred Orientation at the Highest Contrast (degree)');    
    xlim([-dtheta,180+dtheta]);
    for i=lgnmin_e:lgnmax_e
        pick= nP(contrastLevel).ei>0.5 & nLGN==i & nP(contrastLevel).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
        if ~isempty(pick)
            shifts(i-lgnmin_e+1) = sum(abs(aiprefTheta(p2,pick)-aiprefTheta(p1,pick)))/sum(pick)*dtheta;
        end
    end
end


subplot(1,2,2);
bar((lgnmin_e+0.5):(lgnmax_e+0.5),shifts);
xlabel('# LGN input');
ylabel('mean shift in preferred orientation (^{o})');
title('Exc only');
xlim([lgnmin_e,inf]);
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(format,'fig')
        saveas(hOriCmp,[outputfdr,'/','OriCmp-',theme,'.',format]);
    else
        print(hOriCmp,[outputfdr,'/','OriCmp-',theme,'.',format],printDriver,dpi);
    end
end
%% width
if tcReady
    hTCwidth = figure;
    p2 = contrastLevel; p1 = contrastLevel-2;
    pick = nP(contrastLevel).ei>0.5 & nP(p1).pkrate>thres;
    if (0)
    subplot(2,4,1); hold on;
        edges = 0:5:90;
        pickedWidth = width(p1,pick);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','r');
        title('25%');
        xlabel('half width (\theta)')
        ylabel('% exc neurons')
        xlim([0,90]);
    subplot(2,4,2); hold on;
        edges = 0:5:90;
        pickedWidth = width(p2,pick);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','r');
        title('100%');
        xlabel('half width (\theta)')
        ylabel('% exc neurons')
        xlim([0,90]);
    end
    hExcWidth = subplot(1,2,1);
        ctrs = cell(2,1);
        if sum(pick) > 0
            pair1 = width(p1,pick);
            pair2 = width(p2,pick);
            %dTick = 30;
            %tickLim = [0,90];
            %[tick,n0] = autoAxis(tickLim(1),tickLim(2),round(tickLim(2)/dTick),tickLim);
            idTick = 6;
            tick = [0,30,60,90];
            n0 = length(tick);
            tickLabel = num2str(tick');
            lctrsx = (n0-1)*idTick+1;
            lctrsy = lctrsx;
            dctr = (tick(n0)-tick(1))/lctrsx-1;            
            ctrs{1} = linspace(tick(1),tick(n0),lctrsx)+dctr/2;
            ctrs{2} = ctrs{1};
            tickPosY = 0.5:idTick:(lctrsy-1+0.5);
            tickPosX = 0.5:idTick:(lctrsx-1+0.5);
            
            dataPair = [pair1', pair2'];
            denPair = hist3(dataPair,ctrs);
            maxDen = max(max(denPair));
            denPair = denPair/maxDen;
            denPair = denPair(1:(lctrsx-1),1:(lctrsy-1));
            imagesc([1,lctrsx-1],[lctrsy-1,1],denPair');
            hold on
            plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
            title([num2str(sum(pick)/sum(nP(contrastLevel).ei>0.5)*100,'%.1f'),'% qualified neurons']);
            xlabel('half width(25%)')
            ylabel('half width(100%)')
            daspect([1,1,1]);
            
            set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
            colormap(hExcWidth,redOnly);
        end

    pick = nP(contrastLevel).ei<0.5 & nP(p1).pkrate>thres;
    %pick = nP(contrastLevel).ei<0.5 & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(p1).pkrate>thres & nP(p2).pkrate>thres;
    if (0)
    subplot(2,4,5); hold on;
        edges = 0:10:180;
        pickedWidth = width(p1,pick);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','b');
        title('25%');
        xlabel('half width (\theta)')
        ylabel('% inh neurons')
        xlim([0,180]);
    subplot(2,4,6); hold on;
        edges = 0:45:720;
        pickedWidth = width(p2,pick);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','b');
        title('100%');
        xlabel('half width (\theta)')
        ylabel('% inh neurons')
        xlim([0,720]);
    end
    pick = pick & width(p2,:)' < 200 & width(p1,:)' < 200;

    hInhWidth = subplot(1,2,2); 
        idTick = 6;
        ctrs = cell(2,1);
        dTick = 60;
        tickLim = [0,200];

        if sum(pick) > 0
            pair1 = width(p1,pick);
            pair2 = width(p2,pick);
            [tick,n0] = autoAxis(tickLim(1),tickLim(2),round(tickLim(2)/dTick),tickLim);
            tickLabel = num2str(tick');
            lctrsx = (n0-1)*idTick+1;
            lctrsy = lctrsx;
            dctr = (tick(n0)-tick(1))/lctrsx-1;
            ctrs{1} = linspace(tick(1),tick(end),lctrsx)+dctr/2;
            ctrs{2} = ctrs{1};
            tickPosY = 0.5:idTick:(lctrsy-1+0.5);
            tickPosX = 0.5:idTick:(lctrsx-1+0.5);
            
            dataPair = [pair1', pair2'];
            denPair = hist3(dataPair,ctrs);
            maxDen = max(max(denPair));
            denPair = denPair/maxDen;
            denPair = denPair(1:(lctrsx-1),1:(lctrsy-1));
            imagesc([1,lctrsx-1],[lctrsy-1,1],denPair');
            hold on
            plot(lctrsx-0.5:-1:-0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
            title([num2str(sum(pick)/sum(nP(contrastLevel).ei<0.5)*100,'%.1f'),'% qualified neurons']);
            xlabel('half width(25%)')
            ylabel('half width(100%)')
            daspect([1,1,1]);
            
            set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
            colormap(hInhWidth,blueOnly);
        end

    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hTCwidth,[outputfdr,'/','TCwidth-',theme,'.',format]);
        else
            print(hTCwidth,[outputfdr,'/','TCwidth-',theme,'.',format],printDriver,dpi);
        end
    end
end
%% current
    if current
        eoCurr = zeros(contrastLevel,p.nv1e);
        epCurr = zeros(contrastLevel,p.nv1e);
        eoCurrc= zeros(contrastLevel,p.nv1e);
        epCurrc= zeros(contrastLevel,p.nv1e);
        
        ioCurr = zeros(contrastLevel,p.nv1i);
        ipCurr = zeros(contrastLevel,p.nv1i);
        ioCurrc = zeros(contrastLevel,p.nv1i);
        ipCurrc = zeros(contrastLevel,p.nv1i);
        if ~tcReady
            aiprefTheta = zeros(contrastLevel,nv1);
            aiorthTheta = aiprefTheta;
            for i = 1:contrastLevel
                aiprefTheta(i,:) = round(nP(i).prA*180/pi/dtheta)+1;
                bbb = aiprefTheta(i,:)>=ntheta/2+1;
                aiorthTheta(i,bbb) = aiprefTheta(i,bbb)-ntheta/2;
                aiorthTheta(i,~bbb) = aiprefTheta(i,~bbb)+ntheta/2;
            end
        end    

        for j=1:contrastLevel
            nntheta = size(tC(j).currEt,2);
            %pref
            pick = false(ntotal * nntheta,1);
            opPick = aiprefTheta(j,:) + (0:nntheta:((ntotal-1)*nntheta));
            pick(opPick) = true;
            pick = reshape(pick,[nntheta,ntotal])';
            
            eipick = pick;
            eipick(nP(1).ei < 0.5,:) = false;
            epCurr(j,:) = tC(j).currEt(eipick)';
            epCurrc(j,:) = tC(j).currEc(eipick)';
            
            eipick = pick;
            eipick(nP(1).ei > 0.5,:) = false;
            ipCurr(j,:) = tC(j).currEt(eipick)';
            ipCurrc(j,:) = tC(j).currEc(eipick)';
            %orth
            pick = false(ntotal * nntheta,1);
            opPick = aiorthTheta(j,:) + (0:nntheta:((ntotal-1)*nntheta));
            pick(opPick) = true;
            pick = reshape(pick,[nntheta,ntotal])';
        
            eipick = pick;
            eipick(nP(1).ei < 0.5,:) = false;
            eoCurr(j,:) = tC(j).currEt(eipick)';
            eoCurrc(j,:) = tC(j).currEc(eipick)';
        
            eipick = pick;
            eipick(nP(1).ei > 0.5,:) = false;
            ioCurr(j,:) = tC(j).currEt(eipick)';
            ioCurrc(j,:) = tC(j).currEc(eipick)';
        end
        hCurrent = figure;

        for i=1:contrastLevel
            subplot(contrastLevel,2,(i-1)*2+1)
            hold on
            grid on
            errorbar(1,mean(epCurr(i,:)),std(epCurr(i,:)),'Marker','o','Color','r','LineWidth',2,'MarkerSize',7);
            errorbar(2,mean(eoCurr(i,:)),std(eoCurr(i,:)),'Marker','o','Color','r','LineWidth',2,'MarkerSize',7);
            errorbar(3,mean(ipCurr(i,:)),std(ipCurr(i,:)),'Marker','o','Color','b','LineWidth',2,'MarkerSize',7);
            errorbar(4,mean(ioCurr(i,:)),std(ioCurr(i,:)),'Marker','o','Color','b','LineWidth',2,'MarkerSize',7);
            title(['LGN current at ',conLabel(i)]);
            if i==contrastLevel
                set(gca,'XTick',1:4);
                set(gca,'XTickLabel',{'prefE','orthE','prefI','orthI'});
            end
            subplot(contrastLevel,2,(i-1)*2+2)
            hold on
            grid on
            errorbar(1,mean(epCurrc(i,:)),std(epCurrc(i,:)),'Marker','o','Color','r','LineWidth',2,'MarkerSize',7);
            errorbar(2,mean(eoCurrc(i,:)),std(eoCurrc(i,:)),'Marker','o','Color','r','LineWidth',2,'MarkerSize',7);
            errorbar(3,mean(ipCurrc(i,:)),std(ipCurrc(i,:)),'Marker','o','Color','b','LineWidth',2,'MarkerSize',7);
            errorbar(4,mean(ioCurrc(i,:)),std(ioCurrc(i,:)),'Marker','o','Color','b','LineWidth',2,'MarkerSize',7);
            title(['Cortical Exc current at ',conLabel(i)]);
            if i==contrastLevel
                set(gca,'XTick',1:4);
                set(gca,'XTickLabel',{'prefE','orthE','prefI','orthI'});
            end
        end
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(hCurrent,[outputfdr,'/','LGNCurrent-',theme,'.',format]);
            else
                print(hCurrent,[outputfdr,'/','LGNCurrent-',theme,'.',format],printDriver,dpi);
            end
        end
    end
end
    %%    all Tuning curve together
    evenOdd = mod(ntheta,2);
    relTheta = (-(ntheta+evenOdd)/2:(ntheta+evenOdd)/2)*dtheta;
    itheta = ntheta/2 + 1;
    rasterTC = zeros(ntheta+1,ntotal,contrastLevel)-1;
    EmeanTC = zeros(ntheta+1,contrastLevel)-1;
    EquaLTC = zeros(ntheta+1,contrastLevel)-1;
    EquaUTC = zeros(ntheta+1,contrastLevel)-1;
    EmediTC = zeros(ntheta+1,contrastLevel)-1;
    ImeanTC = zeros(ntheta+1,contrastLevel)-1;
    IquaLTC = zeros(ntheta+1,contrastLevel)-1;
    IquaUTC = zeros(ntheta+1,contrastLevel)-1;
    ImediTC = zeros(ntheta+1,contrastLevel)-1;
    TCgLGN = rasterTC;
    TCgLGN_F1 = rasterTC;
    TCgE = rasterTC;
    TCgI = rasterTC;
    EmeanTCgLGN = EmeanTC;
    ImeanTCgLGN = ImeanTC;
    EmeanTCgLGN_F1 = EmeanTC;
    ImeanTCgLGN_F1 = ImeanTC;
    EmeanTCgE = EmeanTC;
    ImeanTCgE = ImeanTC;
    EmeanTCgI = EmeanTC;
    ImeanTCgI = ImeanTC;
    io = {'i','o'};
    for ij = 1:2
        if ~operiod && ij == 2
            continue
        end
        hTCall = figure;
        for i = 1:contrastLevel
            epick = nP(contrastLevel).ei>0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
            ipick = nP(contrastLevel).ei<0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
            for k=1:ntotal    
                %[~,ipA] = max(tC(i).gLGN(k,:).*(1+tC(i).gLGNsc(k,:)));
                % output angle for inh, input for exc
                if ij == 1
                    ipA = round(nP(i).priA(k)*180/pi/dtheta)+1;
                else
                    ipA = round(nP(i).prA(k)*180/pi/dtheta)+1;
                end
                %if k<=p.nv1e
                %    ipA = round(nP(i).priA(k)*180/pi/dtheta)+1;
                %else
                %    ipA = round(nP(i).prA(k)*180/pi/dtheta)+1;
                %end
                if ipA > ntheta
                    ipA = ipA-ntheta;
                end
                %ipA = round(nP(i).prA(k)*180/pi/dtheta)+1;
                if ipA < ntheta/2+1
                    half = [(ipA+ntheta/2):ntheta,1:ipA+ntheta/2];
                else 
                    half = [ipA-ntheta/2:ntheta,1:(ipA-ntheta/2)];
                end
                rho = [tC(i).frate(k,:),tC(i).frate(k,1)];
                gLGN_F1 = [tC(i).gLGN(k,:).*tC(i).gLGNsc(k,:), tC(i).gLGN(k,1)*tC(i).gLGNsc(k,1)];
                gLGN = [tC(i).gLGN(k,:), tC(i).gLGN(k,1)];
                gE = [tC(i).gE(k,:),tC(i).gE(k,1)];
                gI = [tC(i).gI(k,:),tC(i).gI(k,1)];
                rasterTC(:,k,i) = rho(half);
                TCgLGN(:,k,i) = gLGN(half);
                TCgLGN_F1(:,k,i) = gLGN_F1(half);
                TCgE(:,k,i) = gE(half);
                TCgI(:,k,i) = gI(half);
            end
            subplot(contrastLevel,6,(i-1)*6+2)
            hold on
            %raster points
            %plot(relTheta,rasterTC(:,epick,i),'.r','MarkerSize',1);
            % mean TC
            EmeanTC(:,i) = mean(rasterTC(:,epick,i),2);
            EmediTC(:,i) = median(rasterTC(:,epick,i),2);
            EquaLTC(:,i) = EmediTC(:,i)-quantile(rasterTC(:,epick,i),0.25,2);
            EquaUTC(:,i) = quantile(rasterTC(:,epick,i),0.75,2)-EmediTC(:,i);
            EmeanTCgLGN(:,i) = mean(TCgLGN(:,epick,i),2);
            EmeanTCgLGN_F1(:,i) = mean(TCgLGN_F1(:,epick,i),2);
            EmeanTCgE(:,i) = mean(TCgE(:,epick,i),2);
            EmeanTCgI(:,i) = mean(TCgI(:,epick,i),2);
    
            ax1 = gca(gcf);
            axesPos = get(ax1,'Position');
            %ax2 = axes('Position',axesPos,'YAxisLocation','right','Color','None','Box','off','XTick',[]);
            ax2 = axes('Position',axesPos,'YAxisLocation','right','Box','off','XTick',[]);
            axes(ax1);
            hold(ax1,'on');
            errorbar(relTheta, EmediTC(:,i),EquaLTC(:,i),EquaUTC(:,i),'Color','k');
            plot(relTheta,EmeanTC(:,i),'x','Color','k');
            %errorbar(relTheta, EmeanTC(:,i),std(rasterTC(:,epick,i),1,2),'Color','k');
    
            axes(ax2);
            hold(ax2,'on');
            errorbar(relTheta, EmeanTCgLGN(:,i),std(TCgLGN(:,epick,i),1,2),'Color','g');
            errorbar(relTheta, EmeanTCgLGN_F1(:,i),std(TCgLGN_F1(:,epick,i),1,2),'Color','m');
            errorbar(relTheta, EmeanTCgE(:,i),std(TCgE(:,epick,i),1,2),'Color','r');
    
            set(ax1,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'Box','off');
            set(ax2,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'YColor','b');
            axes(ax1);
            set(ax1,'Color','None');
            xticks = -90:30:90;
            xlim([relTheta(1),relTheta(end)]);
            set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
            set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
    
            if i==contrastLevel,xlabel(ax1,'relative angle to pref(^{o})');end
            title([num2str(sum(epick)/p.nv1e*100,'%2.1f'),'% Exc evoked'])
    
            subplot(contrastLevel,6,(i-1)*6+1)
            %errorbar(relTheta, EmeanTCgI(:,i),std(TCgI(:,epick,i),1,2),'Color','b');
            eepick = 1:p.nv1e;
            if i==contrastLevel
                objectMax = EmeanTCgI(itheta,i);
                hold on
                for j=1:contrastLevel-1
                    MeanNormTCgI = TCgI(:,eepick,j)./(ones(size(TCgI,1),1)*TCgI(itheta,eepick,j));
                    MeanNormTCgI = mean(MeanNormTCgI,2);
                    plot(relTheta, MeanNormTCgI*objectMax,LineScheme{j},'Color','b');
                end
                plot(relTheta, EmeanTCgI(:,i),LineScheme{i},'Color','b');
                legend(conLabel);
                title(['Normalized to C',conLabel{i}]);
            else
                errorbar(relTheta, EmeanTCgI(:,i),std(TCgI(:,epick,i),1,2),'Color','b');
            end
            if i==1, legend('gInh');end
            ylabel(['C ',conLabel{i}]);
            xticks = -90:30:90;
            xlim([relTheta(1),relTheta(end)]);
            ylim([0,inf]);
            set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
            set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
    		
    
            subplot(contrastLevel,6,(i-1)*6+5)
            %hold on
            %raster points
            %plot(relTheta,rasterTC(:,ipick,i),'.b','MarkerSize',1);
            % mean TC
            ImeanTC(:,i) = mean(rasterTC(:,ipick,i),2);
            ImediTC(:,i) = median(rasterTC(:,ipick,i),2);
            IquaLTC(:,i) = ImediTC(:,i)-quantile(rasterTC(:,ipick,i),0.25,2);
            IquaUTC(:,i) = quantile(rasterTC(:,ipick,i),0.75,2) - ImediTC(:,i);
            ImeanTCgLGN(:,i) = mean(TCgLGN(:,ipick,i),2);
            ImeanTCgLGN_F1(:,i) = mean(TCgLGN_F1(:,ipick,i),2);
            ImeanTCgE(:,i) = mean(TCgE(:,ipick,i),2);
            ImeanTCgI(:,i) = mean(TCgI(:,ipick,i),2);
    
            ax1 = gca(gcf);
            axesPos = get(ax1,'Position');
            %ax2 = axes('Position',axesPos,'YAxisLocation','right','Color','None','Box','off','XTick',[]);
            ax2 = axes('Position',axesPos,'YAxisLocation','right','Box','off','XTick',[]);
            axes(ax1);
            hold(ax1,'on');
            errorbar(relTheta, ImediTC(:,i),IquaLTC(:,i),IquaUTC(:,i),'Color','k');
            plot(relTheta,ImeanTC(:,i),'x','Color','k');
            %errorbar(relTheta, ImeanTC(:,i),std(rasterTC(:,ipick,i),1,2),'Color','k');
    
            axes(ax2);
            hold(ax2,'on');
            errorbar(relTheta, ImeanTCgLGN(:,i),std(TCgLGN(:,ipick,i),1,2),'Color','g');
            errorbar(relTheta, ImeanTCgLGN_F1(:,i),std(TCgLGN_F1(:,ipick,i),1,2),'Color','m');
            errorbar(relTheta, ImeanTCgE(:,i),std(TCgE(:,ipick,i),1,2),'Color','r');
    
            set(ax1,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'Box','off');
            set(ax2,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'YColor','b');
            axes(ax1);
            set(ax1,'Color','None');
            xticks = -90:30:90;
            xlim([relTheta(1),relTheta(end)]);
            set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
            set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
    
            if i==contrastLevel,xlabel(ax1,'relative angle to pref(^{o})');end
            title([num2str(sum(ipick)/p.nv1i*100,'%2.1f'),'% Inh evoked'])
    
            subplot(contrastLevel,6,(i-1)*6+4)
            errorbar(relTheta, ImeanTCgI(:,i),std(TCgI(:,ipick,i),1,2),'Color','b');
            if i==contrastLevel
                objectMax = ImeanTCgI(itheta,i);
                hold on
                iipick = p.nv1e + (1:p.nv1i);
                for j=1:contrastLevel-1
                    MeanNormTCgI = TCgI(:,iipick,j)./(ones(size(TCgI,1),1)*TCgI(itheta,iipick,j));
                    MeanNormTCgI = mean(MeanNormTCgI,2);
                    plot(relTheta, MeanNormTCgI*objectMax,LineScheme{j},'Color','b');
                end
                plot(relTheta, ImeanTCgI(:,i),'Color','b');
            end
            if i==1, legend('gInh');end
            xticks = -90:30:90;
            xlim([relTheta(1),relTheta(end)]);
            set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
            set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
        end
        
        subplot(2,6,3)
        hold on

        eepick = 1:p.nv1e;
        for i = 1:contrastLevel
            MeanNormTC = rasterTC(:,eepick,i)./(ones(size(rasterTC,1),1)*rasterTC(itheta,eepick,i));
            pick = rasterTC(itheta,eepick,i)>0;
            MeanNormTC = mean(MeanNormTC(:,pick),2);

            mTmp = ones(size(TCgLGN,1),1)*TCgLGN(itheta,eepick,i);
            MeanNormTCgLGN = TCgLGN(:,eepick,i)./mTmp;
            MeanNormTCgLGN = mean(MeanNormTCgLGN,2);

            MeanNormTCgLGN_F1 = TCgLGN_F1(:,eepick,i)./mTmp;
            MeanNormTCgLGN_F1 = mean(MeanNormTCgLGN_F1,2);

            MeanNormTCgE = TCgE(:,eepick,i)./(ones(size(TCgE,1),1)*TCgE(itheta,eepick,i));
            MeanNormTCgE = mean(MeanNormTCgE,2);

            plot(relTheta, MeanNormTCgLGN,LineScheme{i},'Color','g');
            plot(relTheta, MeanNormTCgLGN_F1,LineScheme{i},'Color','m');
            plot(relTheta, MeanNormTCgE,LineScheme{i},'Color','r');
            plot(relTheta, MeanNormTC,LineScheme{i},'Color','k');
        end
        ylabel('normalized to max');
        %yy = ylim();
        %ylim([0,yy(2)]);
        ylim([0,1]);
        xticks = -90:30:90;
        xlim([relTheta(1),relTheta(end)]);
        set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'));
        title('Exc');
    
        subplot(2,6,9)
        ax1 = gca(gcf);
        axesPos = get(ax1,'Position');
        %ax2 = axes('Position',axesPos,'YAxisLocation','right','Color','None','YColor','b');
        ax2 = axes('Position',axesPos,'YAxisLocation','right','YColor','b');
        set(ax1,'Color','None');
        for i=1:contrastLevel
            axes(ax1);
            hold(ax1,'on')
            plot(relTheta, EmeanTC(:,i),LineScheme{i},'Color','k');
            axes(ax2);
            hold(ax2,'on');
            if i == contrastLevel
                h(1) = plot(relTheta, EmeanTCgLGN(:,i),LineScheme{i},'Color','g');
                h(2) = plot(relTheta, EmeanTCgLGN_F1(:,i),LineScheme{i},'Color','m');
                h(3) = plot(relTheta, EmeanTCgE(:,i),LineScheme{i},'Color','r');
            else
                plot(relTheta, EmeanTCgLGN(:,i),LineScheme{i},'Color','g');
                plot(relTheta, EmeanTCgLGN_F1(:,i),LineScheme{i},'Color','m');
                plot(relTheta, EmeanTCgE(:,i),LineScheme{i},'Color','r');
            end
            %plot(relTheta, EmeanTCgI(:,i),LineScheme{i},'Color','b');
        end
        xlabel(ax1,'relative angle to pref(^{o})');
        legend(ax2,h,{'g_{LGN}','gLGN-F_1','g_E'},'Location','NorthEast');
        ylabel(ax1,'Hz'); ylabel(ax2,'Conductance');
        yy = ylim(ax1);
        ylim(ax1,[0,1.3*yy(2)]);
        yy = ylim(ax2);
        ylim(ax2,[0,1.3*yy(2)]);
        axes(ax1);
        set(ax1,'Color','None');
        xticks = -90:30:90;
        xlim(ax1,[relTheta(1),relTheta(end)]);
        xlim(ax2,[relTheta(1),relTheta(end)]);
        set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
        set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
    
        subplot(2,6,6)
        hold on
        
        iipick = p.nv1e + (1:p.nv1i);
        for i = 1:contrastLevel
            MeanNormTC = rasterTC(:,iipick,i)./(ones(size(rasterTC,1),1)*rasterTC(itheta,iipick,i));
            pick = rasterTC(itheta,iipick,i)>0;
            MeanNormTC = mean(MeanNormTC(:,pick),2);

            mTmp = ones(size(TCgLGN,1),1)*TCgLGN(itheta,iipick,i);
            MeanNormTCgLGN = TCgLGN(:,iipick,i)./mTmp;
            MeanNormTCgLGN = mean(MeanNormTCgLGN,2);

            MeanNormTCgLGN_F1 = TCgLGN_F1(:,iipick,i)./mTmp;
            MeanNormTCgLGN_F1 = mean(MeanNormTCgLGN_F1,2);

            MeanNormTCgE = TCgE(:,iipick,i)./(ones(size(TCgE,1),1)*TCgE(itheta,iipick,i));
            MeanNormTCgE = mean(MeanNormTCgE,2);

            plot(relTheta, MeanNormTCgLGN,LineScheme{i},'Color','g');
            plot(relTheta, MeanNormTCgLGN_F1,LineScheme{i},'Color','m');
            plot(relTheta, MeanNormTCgE,LineScheme{i},'Color','r');
            plot(relTheta, MeanNormTC,LineScheme{i},'Color','k');
        end
        ylabel('normalized to max');
        %yy = ylim();
        %ylim([0,yy(2)]);
        ylim([0,1]);
        title('Inh');
        xticks = -90:30:90;
        xlim([relTheta(1),relTheta(end)]);
        set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'));

        subplot(2,6,12)
        ax1 = gca(gcf);
        axesPos = get(ax1,'Position');
        %ax2 = axes('Position',axesPos,'YAxisLocation','right','Color','None');
        ax2 = axes('Position',axesPos,'YAxisLocation','right');
        axes(ax1);
        for i=1:contrastLevel
            hold(ax1,'on')
            plot(relTheta, ImeanTC(:,i),LineScheme{i},'Color','k');
        end
        axes(ax2);
        hold(ax2,'on');
        for i=1:contrastLevel
            plot(relTheta, ImeanTCgLGN(:,i),LineScheme{i},'Color','g');
            plot(relTheta, ImeanTCgLGN_F1(:,i),LineScheme{i},'Color','m');
            plot(relTheta, ImeanTCgE(:,i),LineScheme{i},'Color','r');
            %plot(relTheta, EmeanTCgI(:,i),LineScheme{i},'Color','b');
        end
        set(ax2,'YColor','b')
        xlabel(ax1,'relative angle to pref(^{o})');
        ylabel(ax1,'Hz'); ylabel(ax2,'Conductance');
        yy = ylim(ax1);
        ylim(ax1,[0,yy(2)]);
        yy = ylim(ax2);
        ylim(ax2,[0,yy(2)]);
        axes(ax1);
        set(ax1,'Color','None');
        xticks = -90:30:90;
        xlim(ax1,[relTheta(1),relTheta(end)]);
        xlim(ax2,[relTheta(1),relTheta(end)]);
        set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
        set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
        if ~isempty(format)
            set(gcf,'Renderer','Painters')
            %set(gcf,'Renderer','zbuffer')
            %set(gcf,'Renderer','OpenGL')
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
            if strcmp(format,'fig')
                saveas(hTCall,[outputfdr,'/','TCall',io{ij},'-',theme,'.',format]);
            else
                print(hTCall,[outputfdr,'/','TCall',io{ij},'-',theme,'.',format],printDriver,dpi);
            end
        end
        %% TC of types
        hTypeTC = figure;

        ntypes = p.ntypeE+p.ntypeI;
        meanTC = zeros(ntheta+1,contrastLevel,ntypes);
        meanTCgLGN = meanTC;
        meanTCgE = meanTC;
        for i = 1:contrastLevel
            epick = nP(contrastLevel).ei>0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
            ipick = nP(contrastLevel).ei<0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
            for j = 1:p.ntypeE
                pick = epick & [p.typeE == j; false(p.nv1i,1)];
                subplot(contrastLevel+1,ntypes,(i-1)*ntypes+j)
                hold on 
                meanTC(:,i,j) = mean(rasterTC(:,pick,i),2);
                meanTCgLGN(:,i,j) = mean(TCgLGN(:,pick,i),2);
                meanTCgLGN_F1(:,i,j) = mean(TCgLGN_F1(:,pick,i),2);
                meanTCgE(:,i,j) = mean(TCgE(:,pick,i),2);
                %[ax, h1, h2] = plotyy(meanTC(:,i,j),std(rasterTC(:,pick,i),1,2),...
                %                [meanTCgLGN(:,i,j),meanTCgE(:,i,j)],[std(TCgLGN(:,pick,i),1,2),std(TCgE(:,pick,i),1,2)],...
                %                'errorbar');
                ax1 = gca(gcf);
                axesPos = get(ax1,'Position');
                ax2 = axes('Position',axesPos,'YAxisLocation','right','Color','None','Box','off','XTick',[]);
                axes(ax1);
                errorbar(relTheta, meanTC(:,i,j),std(rasterTC(:,pick,i),1,2),'Color','k');

                axes(ax2);
                hold(ax2,'on');
                errorbar(relTheta, meanTCgLGN(:,i,j),std(TCgLGN(:,pick,i),1,2),'Color','g');
                errorbar(relTheta, meanTCgLGN_F1(:,i,j),std(TCgLGN_F1(:,pick,i),1,2),'Color','m');
                errorbar(relTheta, meanTCgE(:,i,j),std(TCgE(:,pick,i),1,2),'Color','r');

                set(ax1,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'Box','off');
                set(ax2,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'YColor','b');

                evoked = [num2str(100*sum(pick)/sum(p.typeE == j),'%3.1f'),'%'];
                if i==contrastLevel && j==1,ylabel(ax1,'Hz');ylabel(ax2,'Conductance');end
                if i==1 
                    title([evoked,' ' ,p.Etypes{j}]);
                else
                    title(evoked);
                end

            end
            for j = 1:p.ntypeI
                pick = ipick & [false(p.nv1e,1); p.typeI == j];
                subplot(contrastLevel+1,ntypes,(i-1)*ntypes+p.ntypeE+j)
                hold on
                meanTC(:,i,p.ntypeE+j) = mean(rasterTC(:,pick,i),2);
                meanTCgLGN(:,i,p.ntypeE+j) = mean(TCgLGN(:,pick,i),2);
                meanTCgLGN_F1(:,i,p.ntypeE+j) = mean(TCgLGN_F1(:,pick,i),2);
                meanTCgE(:,i,p.ntypeE+j) = mean(TCgE(:,pick,i),2);

                ax1 = gca(gcf);
                axesPos = get(ax1,'Position');
                ax2 = axes('Position',axesPos,'YAxisLocation','right','Color','None','Box','off','XTick',[]);
                axes(ax1);
                errorbar(relTheta, meanTC(:,i,p.ntypeE+j),std(rasterTC(:,pick,i),1,2),'Color','k');

                axes(ax2);
                hold(ax2,'on');
                errorbar(relTheta, meanTCgLGN(:,i,p.ntypeE+j),std(TCgLGN(:,pick,i),1,2),'Color','g');
                errorbar(relTheta, meanTCgLGN_F1(:,i,p.ntypeE+j),std(TCgLGN_F1(:,pick,i),1,2),'Color','m');
                errorbar(relTheta, meanTCgE(:,i,p.ntypeE+j),std(TCgE(:,pick,i),1,2),'Color','r');

                set(ax1,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'Box','off');
                set(ax2,'XLim',[relTheta(1),relTheta(end)],'YLim',[0,inf],'YColor','b');

                evoked = [num2str(100*sum(pick)/sum(p.typeI == j),'%3.1f'),'%'];
                if i==1 
                    title([evoked,' ' ,p.Itypes{j}]);
                else
                    title(evoked);
                end
            end
        end
        for j=1:ntypes
            subplot(contrastLevel+1,ntypes,contrastLevel*ntypes+j)
            hold on
            maxMeanTC = ones(ntheta+1,1) * max(meanTC(:,:,j));
            maxMeanTCgLGN = ones(ntheta+1,1) * max(meanTCgLGN(:,:,j));
            maxMeanTCgLGN_F1 = ones(ntheta+1,1) * max(meanTCgLGN_F1(:,:,j));
            maxMeanTCgE = ones(ntheta+1,1) * max(meanTCgE(:,:,j));
            for i = 1:contrastLevel
                plot(relTheta,meanTC(:,i,j)./maxMeanTC(:,i),LineScheme{i},'Color','k');
                plot(relTheta,meanTCgLGN(:,i,j)./maxMeanTCgLGN(:,i),LineScheme{i},'Color','g');
                plot(relTheta,meanTCgLGN_F1(:,i,j)./maxMeanTCgLGN_F1(:,i),LineScheme{i},'Color','m');
                plot(relTheta,meanTCgE(:,i,j)./maxMeanTCgE(:,i),LineScheme{i},'Color','r');
            end
            ylabel('normalized to max');
            xlabel('relative angle(^{o})');
            ylim([0,inf]);
            if j==1
                title('All contrasts');
            end    
        end

        if ~isempty(format)
            set(gcf,'Renderer','Painters')
            %set(gcf,'Renderer','zbuffer')
            %set(gcf,'Renderer','OpenGL')
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition.*2.5);
            if strcmp(format,'fig')
                saveas(hTypeTC,[outputfdr,'/','TypeTC',io{ij},'-',theme,'.',format]);
            else
                print(hTypeTC,[outputfdr,'/','TypeTC',io{ij},'-',theme,'.',format],printDriver,dpi);
            end
        end
    end

    hFRThetaHeat = figure;
    ORF = p.typeE==1;
    SRF = p.typeE==2;
    Inh = p.nv1e+(1:p.nv1i);

    dddtheta = pi/ntheta;
    nfr = 3;
    ctrs = cell(2,1);
    ctrs{1} = 0:dddtheta:pi;
    lctrsx = length(ctrs{1});
    lctrsy = nfr+1;
    dTickX = 3/ntheta;
    tickPosX = 0.5:lctrsx*dTickX:lctrsx+0.5;
    dTickY = 0.2;
    %tickPosY = 0.5:lctrsy*dTickY:lctrsy+0.5;
    tickLabelX = num2str((linspace(0,pi,length(tickPosX)))'.*180/pi);
    for i = 1:contrastLevel

        hExc = subplot(contrastLevel,3,(i-1)*3+1);
        frTarget = nP(i).pkrate(ORF);
        thetaTarget = nP(i).prA(ORF);
        maxfr = ceil(max(frTarget));
        minfr = floor(min(frTarget));
        if minfr~=maxfr
        ctrs{2} = linspace(minfr,maxfr,lctrsy);
        [tick,n0] = autoAxis(minfr,maxfr,round(1/dTickY),[0,inf]);
        tickLabelY = num2str(flipud(tick'));
        lctrsy = (n0-1)*nfr+1;
        ctrs{2} = linspace(tick(1),tick(n0),lctrsy);
        tickPosY = 0.5:nfr:(lctrsy-1+0.5);
        %tickLabelY = flipud(num2str((linspace(minfr,maxfr,length(tickPosY)))'));
        pair = [thetaTarget,frTarget];
        den = hist3(pair,ctrs);
        den = den(1:ntheta,1:lctrsy-1);
        mdeny = max(den);
        [mden, ix] = max(mdeny);
        den = den/mden;
        imagesc([1,lctrsx-1],[lctrsy-1,1],den');
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        title(['ORF maxFR ', num2str(max(frTarget)),' \theta = ',num2str(ix*dddtheta*180/pi)]);
        xlabel('\theta');
        ylabel('FR Hz');
        colormap(hExc,redOnly);
        end

        hExc = subplot(contrastLevel,3,(i-1)*3+2);
        frTarget = nP(i).pkrate(SRF);
        thetaTarget = nP(i).prA(SRF);
        maxfr = ceil(max(frTarget));
        minfr = floor(min(frTarget));
        if minfr~=maxfr
        ctrs{2} = linspace(minfr,maxfr,lctrsy);
        [tick,n0] = autoAxis(minfr,maxfr,round(1/dTickY),[0,inf]);
        tickLabelY = num2str(flipud(tick'));
        lctrsy = (n0-1)*nfr+1;
        ctrs{2} = linspace(tick(1),tick(n0),lctrsy);
        tickPosY = 0.5:nfr:(lctrsy-1+0.5);
        %tickLabelY = flipud(num2str((linspace(minfr,maxfr,length(tickPosY)))'));
        pair = [thetaTarget,frTarget];
        den = hist3(pair,ctrs);
        den = den(1:ntheta,1:lctrsy-1);
        mdeny = max(den);
        [mden, ix] = max(mdeny);
        den = den/mden;
        imagesc([1,lctrsx-1],[lctrsy-1,1],den');
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        title(['SRF maxFR ', num2str(max(frTarget)),' \theta = ',num2str(ix*dddtheta*180/pi)]);
        xlabel('\theta');
        ylabel('FR Hz');
        colormap(hExc,redOnly);
        end

        hInh = subplot(contrastLevel,3,(i-1)*3+3);
        frTarget = nP(i).pkrate(Inh);
        thetaTarget = nP(i).prA(Inh);
        maxfr = ceil(max(frTarget));
        minfr = floor(min(frTarget));
        if minfr~=maxfr
        [tick,n0] = autoAxis(minfr,maxfr,round(1/dTickY),[0,inf]);
        lctrsy = (n0-1)*nfr+1;
        tickLabelY = num2str(flipud(tick'));
        tickPosY = 0.5:nfr:(lctrsy-1+0.5);
        ctrs{2} = linspace(tick(1),tick(n0),lctrsy);
        %tickLabelY = flipud(num2str((linspace(minfr,maxfr,length(tickPosY)))'));
        pair = [thetaTarget,frTarget];
        den = hist3(pair,ctrs);
        den = den(1:ntheta,1:lctrsy-1);
        mdeny = max(den);
        [mden, ix] = max(mdeny);
        den = den/mden;
        imagesc([1,lctrsx-1],[lctrsy-1,1],den');
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        title(['Inh maxFR ', num2str(max(frTarget)),' \theta = ',num2str(ix*dddtheta*180/pi)]);
        xlabel('\theta');
        ylabel('FR Hz');
        colormap(hInh,blueOnly);
        end
    end
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hFRThetaHeat,[outputfdr,'/','FRThetaHeat-',theme,'.',format]);
        else
            print(hFRThetaHeat,[outputfdr,'/','FRThetaHeat-',theme,'.',format],printDriver,dpi);
        end
    end
    
disp('population statistics done');
toc;
return
end
function [X, Y, Z] = multiLGNspatialKernel(posx,posy,x,y,n,s,scale)
    Aa = 14.88 * (180/pi)^2;
    Ab = Aa*0.97;
    rc = 5.61; %deg
    rs = 16.98*rc/5.61; %deg
    siga = rc/sqrt(2);
    sigb = rs/sqrt(2);
    Z = zeros(length(y),length(x));
    [X, Y] = meshgrid(x,y);
    for i = 1:n
        Z = Z + s(i)*(Aa/(siga^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/siga^2/2)...
            - Ab/(sigb^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/sigb^2/2));
    end
    Z = Z/scale;
end
