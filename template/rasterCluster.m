function clusterFilename = rasterCluster(theme,lgnfile,eps,ctheta,conMatfile,coMatfile,sfile,clusterFile,nthres,cthres,Mov,format,processed,clusterOnly)
    if ~clusterOnly
        twindow = 5; %s
        t1 = 21; %s
        t0 = t1-twindow;  %s
        timebin = 10 /1000; % ms
        nTestCluster = 10;
        maxTau = 25; % in timebin
    end
    pPosition = [0, 0, 1280, 720];
    mSize = 8;
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r100';
    end

FontSize = 10;
set(0,'DefaultAxesFontSize',FontSize)
    load(lgnfile,'p','nSubLGN','v1Map','LGNpos');
    %load(sfile);
    fid = fopen(conMatfile);
    n = p.nv1e;
    m = fread(fid,[p.nv1,p.nv1],'int8');
    mee = m(1:n,1:n);
    fclose(fid);
     
    load(coMatfile);
    [~,sid] = sort(mee,'descend');
    position = sid(1:nthres,:);
    if ~clusterOnly
        epsStr = num2str(eps,'%04d');
        thetaStr = num2str(ctheta,'%02d');
        if processed
            DIR = [theme,'/spike_wise/',epsStr,'/']; 
            [tspI,l] = readSpikes(DIR,p.nv1,[thetaStr,'-spikes.dat']);
        else
            DIR = [theme,'/',epsStr,'/',thetaStr,'/']; 
            [tspI,l] = readSpikes(DIR,p.nv1,['spikes.dat']);
        end
        disp('spikes read');
    end
    if isempty(clusterFile)
        tic;
%         poolobj =gcp('nocreate');
%         delete(poolobj);
%         parpool(nthres);
        seq = randperm(n);
        clusterID = zeros(n,1);
        cluster = cell(n,1);
        cCore = cell(n,1);
        clusterSize = zeros(n,1);
        nCluster = 1;
        cluster{1} = seq(1);
        clusterSize(1) = 1;
        clusterID(seq(1)) = 1;
        for k=2:n
            i = seq(k);
            joined = false;
            for j=1:nCluster
                if clusterSize(j) == 1
                    overlap = inset(position(:,i),position(:,cluster{j}));
                    if  sum(overlap) >= cthres
                        cCore{j} = position(overlap,i);
                        joined = true;
                    end
                else
                    nOverlap = sum(inset(position(:,i),cCore{j}));
                    if nOverlap >= cthres
                        joined = true;
                    end
                end
                if joined
                    cluster{j} = [cluster{j};i];
                    clusterSize(j) = clusterSize(j) + 1;
                    clusterID(i) = j;
                    break;
                end
            end
            if ~joined
                nCluster = nCluster + 1;
                cluster{nCluster} = i;
                clusterSize(nCluster) = clusterSize(nCluster) + 1;
                clusterID(i) = nCluster;
            end
            if mod(k,ceil(n*0.1)) == 0
                disp([num2str(k/n*100),'%']);
            end
        end
        [clusterSize,sortedIndC] = sort(clusterSize,'descend');
        cluster = cluster(sortedIndC);
        cCore = cCore(sortedIndC);
        ylims = zeros(nCluster,2);
        ylims(:,2) = cumsum(clusterSize(1:nCluster))+0.5;
        ylims(:,1) = cumsum([0;clusterSize(1:(nCluster-1))])+0.5;
        disp('clusterized');
        clusterFilename = ['Cluster-',theme,'-',num2str(nthres),'-',num2str(cthres),'.mat'];
        save(clusterFilename,'cluster','cCore','nCluster','seq','clusterSize','ylims','clusterID','nthres','cthres');
        toc;
    else
        clusterFilename = clusterFile;
        load(clusterFile);
        disp('clusters loaded');
    end
    disp(['mean cluster size = ',num2str(mean(clusterSize(1:nCluster)))]);
    disp([num2str(nCluster),' clusters']);
    clear mee;
    if ~clusterOnly
        id = 0;
        h = figure;
%         subplot(2,1,1);
        hold on
        hPlot = 1;
        for i=1:nCluster
            iCluster = cluster{i};
            nEvents = sum(l(iCluster));
            if nEvents>0
                isp = zeros(nEvents,1);
                tsp = zeros(nEvents,1);
                iEvent = 0;
                for j=1:clusterSize(i)
                    if l(iCluster(j))>0
                        cEvent = 1:l(iCluster(j));
                        isp(iEvent+cEvent) = id + j; 
                        tsp(iEvent+cEvent) = tspI(iCluster(j)).tsp;
                        iEvent = iEvent + l(iCluster(j));
                    end
                end
                eRaster(hPlot) = plot(tsp*1000,isp,'.','MarkerSize',mSize); 
                hPlot = hPlot + 1;
            end
            id = id + clusterSize(i);
        end 
        iseq=n+1:p.nv1;
        nEvents = sum(l(iseq));
        isp = zeros(nEvents,1);
        tsp = zeros(nEvents,1);
        iEvent = 0;
        for i = iseq
            cEvent = 1:l(i);
            isp(iEvent+cEvent) = i;
            tsp(iEvent+cEvent) = tspI(i).tsp;
            iEvent = iEvent + l(i);
        end
        plot(tsp*1000,isp,'.','MarkerSize',mSize);
        xlabel('ms'); ylabel('clustered Neuron ID');
        title([theme,'-',epsStr,'-',thetaStr]);
        %ylim([0,300]);
        xlim([t0,t1]*10^3);
        %% 
%         v1pos = zeros(n,2);
%         for i=1:n
%             for j =1:length(nSubLGN{i})
%                 v1pos(i,:) = v1pos(i,:) + sum(LGNpos(abs(v1Map{i,j}),:),1);
%             end
%             v1pos(i,:) = v1pos(i,:)/sum(abs(nSubLGN{i}));
%         end
% 
%         if Mov
%             h2 = subplot(2,2,3);
%         else
%             subplot(2,2,3);
%             hold on
%             for i=1:nCluster
%                 plot(v1pos(cluster{i},1),v1pos(cluster{i},2),'.','MarkerSize',mSize);
%             end
%             axis equal
        if ~isempty(format)
            figname = ['raster-',theme,'-',epsStr,'-',thetaStr,'.',format];
            if strcmp(format,'fig')
                savefig(h,figname);
            else
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                print(h,['raster-',theme,'-',epsStr,'-',thetaStr,'.',format],printDriver,dpi);
            end
        end
%             %% xcorr
% %             close(h);
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
         title([epsStr,thetaStr]);
 
         if ~isempty(format)
             figname = ['Xcorr-',theme,'-',epsStr,'-',thetaStr,'_',numstr(nthres),'-',num2str(cthres),'.',format];
             if strcmp(format,'fig')
                 savefig(h,figname);
             else
                 set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                 print(h,figname,printDriver,dpi);
             end
         end
% %         end
    end
end
function mark = inset(a,b)
    n = length(a);
    mark = false(n,1);
    for i=1:n;
        if ~isempty(find(b==a(i),1)),mark(i)=true;end
    end
end
