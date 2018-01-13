format = 'png';
pPosition = [0, 0, 1280, 720];
load('Cluster-mE4');
load('1xu-n3-s911.mat');
theme = 'mE4-12-30';
eps =0125;ctheta=2;    
    epsStr = num2str(eps,'%04d');
    thetaStr = num2str(ctheta,'%02d');
    DIR = [theme,'/spike_wise/',epsStr,'/']; 
    [tspI,l] = readSpikes(DIR,p.nv1,[thetaStr,'-spikes.dat']);
t0 = 1;  %s
t1 = 21; %s
timebin = 5 /1000; % ms
% nTestCluster = 10;
maxTau = 100; % in timebin
h=figure;
%             subplot(2,2,4);
    hold on
    nt = nearest((t1-t0)/timebin);

%             iTestCluster = randperm(nCluster,nTestCluster);
    iTestCluster = 1:nCluster;
    nxcluster = 100;
    nTestCluster = length(iTestCluster);
    xcorrVec = zeros(maxTau*2+1,nTestCluster);
    tranges = linspace(t0,t1,nt+1);
    tauRanges = (-(maxTau):(maxTau))*timebin;
    CbinedSpikeTrains = zeros(nt,nTestCluster);
    for i=1:nTestCluster
        j = iTestCluster(i);
        sC = clusterSize(j);
        binedSpikeTrains = zeros(nt,sC);
        if sC > 1
            noneZero = true(sC,1);
            for k = 1:clusterSize(j)
                 if l(cluster{j}(k)) > 0
                     tempCounts = histc(tspI(cluster{j}(k)).tsp,tranges);
                     binedSpikeTrains(:,k) = tempCounts(1:nt);
                 else
                     noneZero(k) = false;
                 end
            end
            if sum(noneZero)>1
                CbinedSpikeTrains(:,i) = mean(binedSpikeTrains,2); 
                xcorrMat = xcorr(binedSpikeTrains(:,noneZero),maxTau,'coeff');
                sC = sC - sum(~noneZero);
                autoColumns = (1:sC);
                autoColumns = (autoColumns-1)*sC +autoColumns;
                column = true(sC*sC,1);
                column(autoColumns) = false;
                xcorrVec(:,i) = mean(xcorrMat(:,column),2);
            end
        end
    end
    xcorrVecMean = mean(xcorrVec,2);
    xcorrVecStd = std(xcorrVec,1,2);
    errorbar(tauRanges,xcorrVecMean,xcorrVecStd,'-','LineWidth',2);
%             plot(tauRanges,xcorrVec,'--');

    [~,selectedClusters] = sort(max(xcorrVec),'descend');
    selectedClusters = selectedClusters(1:nxcluster);
    clear xcorrVec;
    xcorrMat = xcorr(CbinedSpikeTrains(:,selectedClusters),maxTau,'coeff');
    autoColumns = (1:nxcluster);
    autoColumns = (autoColumns-1)*nxcluster +autoColumns;
    column = true(nxcluster*nxcluster,1);
    column(autoColumns) = false;
    xcorrCVecMean = mean(xcorrMat(:,column),2);
    xcorrCVecStd = std(xcorrMat(:,column),1,2);
    errorbar(tauRanges,xcorrCVecMean,xcorrCVecStd,':','LineWidth',2);

    ylim([0,inf]);
    legend({'intraCluster','interCluster'});
    xlabel('\tau s');
    ylabel('xcorr');
    figname = ['Xcorr-',theme,'-',epsStr,'-',thetaStr,'.',format];
    title(figname);
    
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r150';
    end
    if ~isempty(format)
        figname = ['Xcorr-',theme,'-',epsStr,'-',thetaStr,'.',format];
        if strcmp(format,'fig')
            saveas(h,figname);
        else
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(h,figname,printDriver,dpi);
        end
    end
