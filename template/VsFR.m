function VsFR(lgnfile,theme,ntheta,contrastLevel,ndperiod,loadData,format,dir,threads)
    if nargin < 9
        threads = 1;
        if nargin < 8
            dir = theme;
            if nargin < 7
                format ='fig';
                if nargin < 6
                    loadData = false;
                    if nargin < 5
                        ndperiod = 25;
                        if nargin < 4
                            contrastLevel = 4;
                            if nargin < 3
                                ntheta = 12;
                            end
                        end
                    end
                end
            end
        end
    end
	if threads > 1
		pool = gcp('nocreate');
		if ~isempty(pool)
			if pool.NumWorkers ~= threads  
				delete(pool);
				pool = parpool(threads);
			end
		else 
			pool = parpool(threads);
		end
	end
    
    load(lgnfile,'p');
    epick = 1:p.nv1e;
    ipick = (p.nv1e+1):p.nv1;
    ntotal = p.nv1;
    dir = [dir,'/'];
    VsEdgeE = 0:0.1:1.2;
    VsStdEdgeE = 0:0.05:0.5;
    VsScEdgeE = 0:0.1:2;
    VsEdgeI = 0:0.1:2.0;
    VsStdEdgeI = 0:0.1:1.0;
    VsScEdgeI = 0:0.1:2;
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r300';
    end
    maxPeriodi = (ndperiod-1)*1/4+1;
    maxThetai = ntheta/2;
    if ~loadData
        priA = zeros(p.nv1,contrastLevel);
        Vs = zeros(p.nv1,ntheta,contrastLevel);
        VsSTD = Vs;
        VsSC = Vs;
        rate = Vs;
        for i=1:contrastLevel
            disp(['loading c',num2str(i)]);
            eps = 125*2^(i-1);
            conLabel(i) = {[num2str(eps/10,'%.1f'),'%']};
            DIR = [theme,'/',num2str(eps,'%04d')];
            [priA(:,i),Vs(:,:,i),VsSTD(:,:,i),VsSC(:,:,i),rate(:,:,i),aC(i)] = readVsFR(DIR,ntheta,ntotal,ndperiod);
            disp('loaded, recalibrating');
            for j =1:ntotal
                if priA(j,i) < maxThetai+1
                    half = [(priA(j,i)+maxThetai):ntheta,1:(priA(j,i)+maxThetai-1)];
                else 
                    half = [priA(j,i)-maxThetai:ntheta,1:(priA(j,i)-maxThetai-1)];
                end
                Vs(j,:,i) = Vs(j,half,i);
                VsSC(j,:,i) = VsSC(j,half,i);
                VsSTD(j,:,i) = VsSTD(j,half,i);
                rate(j,:,i) = rate(j,half,i);
                aC(i).gtot(j,:,:)    = aC(i).gtot(j,:,half);
                aC(i).Itot(j,:,:)    = aC(i).Itot(j,:,half);
                aC(i).Vs(j,:,:)      = aC(i).Vs(j,:,half);
                aC(i).gLGN(j,:,:)    = aC(i).gLGN(j,:,half);
                aC(i).gE(j,:,:)      = aC(i).gE(j,:,half);
                aC(i).gI(j,:,:)      = aC(i).gI(j,:,half);
                aC(i).gEn(j,:,:)     = aC(i).gEn(j,:,half);
                aC(i).gIn(j,:,:)     = aC(i).gIn(j,:,half);
                aC(i).gP(j,:,:)      = aC(i).gP(j,:,half);
                aC(i).V(j,:,:)       = aC(i).V(j,:,half);
                aC(i).gtotstd(j,:,:) = aC(i).gtotstd(j,:,half);
                aC(i).Itotstd(j,:,:) = aC(i).Itotstd(j,:,half);
                aC(i).VsSTD(j,:,:)   = aC(i).VsSTD(j,:,half);
                aC(i).gLGNstd(j,:,:) = aC(i).gLGNstd(j,:,half);
                aC(i).gEstd(j,:,:)   = aC(i).gEstd(j,:,half);
                aC(i).gIstd(j,:,:)   = aC(i).gIstd(j,:,half);
                aC(i).gEnstd(j,:,:)  = aC(i).gEnstd(j,:,half);
                aC(i).gInstd(j,:,:)  = aC(i).gInstd(j,:,half);
                aC(i).gPstd(j,:,:)   = aC(i).gPstd(j,:,half);
                aC(i).Vstd(j,:,:)    = aC(i).Vstd(j,:,half);
                for k=1:ntheta
                   [~,maxIperiod] = max(aC(i).Vs(j,:,k));
                    if maxIperiod <= maxPeriodi
                        iperiod = [(maxIperiod+(ndperiod-1)*3/4):ndperiod,1:(maxIperiod+(ndperiod-1)*3/4-1)];
                    else
                        iperiod = [(maxIperiod-maxPeriodi):ndperiod,1:(maxIperiod-maxPeriodi-1)];
                    end
                    aC(i).gtot(j,:,k)    = aC(i).gtot(j,iperiod,k);
                    aC(i).Itot(j,:,k)    = aC(i).Itot(j,iperiod,k);
                    aC(i).Vs(j,:,k)      = aC(i).Vs(j,iperiod,k);
                    aC(i).gLGN(j,:,k)    = aC(i).gLGN(j,iperiod,k);
                    aC(i).gE(j,:,k)      = aC(i).gE(j,iperiod,k);
                    aC(i).gI(j,:,k)      = aC(i).gI(j,iperiod,k);
                    aC(i).gEn(j,:,k)     = aC(i).gEn(j,iperiod,k);
                    aC(i).gIn(j,:,k)     = aC(i).gIn(j,iperiod,k);
                    aC(i).gP(j,:,k)      = aC(i).gP(j,iperiod,k);
                    aC(i).V(j,:,k)       = aC(i).V(j,iperiod,k);
                    aC(i).gtotstd(j,:,k) = aC(i).gtotstd(j,iperiod,k);
                    aC(i).Itotstd(j,:,k) = aC(i).Itotstd(j,iperiod,k);
                    aC(i).VsSTD(j,:,k)   = aC(i).VsSTD(j,iperiod,k);
                    aC(i).gLGNstd(j,:,k) = aC(i).gLGNstd(j,iperiod,k);
                    aC(i).gEstd(j,:,k)   = aC(i).gEstd(j,iperiod,k);
                    aC(i).gIstd(j,:,k)   = aC(i).gIstd(j,iperiod,k);
                    aC(i).gEnstd(j,:,k)  = aC(i).gEnstd(j,iperiod,k);
                    aC(i).gInstd(j,:,k)  = aC(i).gInstd(j,iperiod,k);
                    aC(i).gPstd(j,:,k)   = aC(i).gPstd(j,iperiod,k);
                    aC(i).Vstd(j,:,k)    = aC(i).Vstd(j,iperiod,k);
                end       
            end
            disp('finished');
        end
        save([theme,'-VsFR.mat'],'priA','Vs','VsSTD','VsSC','rate','aC','-v7.3');
    else 
        load([theme,'-VsFR.mat'],'priA','Vs','VsSTD','VsSC','rate','aC');
    end
    % Vs to Firing Rate transfer function averaged over period
    for k=1:ntheta
        hAvPeriodHeatVsVsSTD = figure;
        hAvPeriodHeatVsVsSC = figure;
        hVsDist = figure;
        k
        for i=1:contrastLevel
            i
            figure(hAvPeriodHeatVsVsSTD);
            ax1 = subplot(4,4,(i-1)*4+1);
            ax2 = subplot(4,4,(i-1)*4+2);
                [binRateE,binRateStdE,VsID] = binVsFr(Vs(epick,k,i),VsSTD(epick,k,i),rate(epick,k,i),VsEdgeE,VsStdEdgeE);
                heatVsStdFr(VsEdgeE,VsStdEdgeE,binRateE,binRateStdE,ax1,ax2,'VsSTD');
            figure(hAvPeriodHeatVsVsSC);
            ax3 = subplot(4,4,(i-1)*4+1);
            ax4 = subplot(4,4,(i-1)*4+2);
                [binRateE,binRateStdE,~] = binVsFr(Vs(epick,k,i),VsSC(epick,k,i),rate(epick,k,i),VsEdgeE,VsScEdgeE,VsID);
                heatVsStdFr(VsEdgeE,VsScEdgeE,binRateE,binRateStdE,ax3,ax4,'VsSC');

            figure(hAvPeriodHeatVsVsSTD);
            ax1 = subplot(4,4,(i-1)*4+3);
            ax2 = subplot(4,4,(i-1)*4+4);
                [binRateI,binRateStdI,VsID] = binVsFr(Vs(ipick,k,i),VsSTD(ipick,k,i),rate(ipick,k,i),VsEdgeI,VsStdEdgeI);
                heatVsStdFr(VsEdgeI,VsStdEdgeI,binRateI,binRateStdI,ax1,ax2,'VsSTD');
            figure(hAvPeriodHeatVsVsSC);
            ax3 = subplot(4,4,(i-1)*4+3);
            ax4 = subplot(4,4,(i-1)*4+4);
                [binRateI,binRateStdI,~] = binVsFr(Vs(ipick,k,i),VsSC(ipick,k,i),rate(ipick,k,i),VsEdgeI,VsScEdgeI,VsID);
                heatVsStdFr(VsEdgeI,VsScEdgeI,binRateI,binRateStdI,ax3,ax4,'VsSC');
        end
        if ~isempty(format)
            fname = [theme,'/VsFRVsSCheat-',theme,'_',num2str(k-maxThetai,'%02d'),'.',format];
            if strcmp(format,'fig')
                saveas(hAvPeriodHeatVsVsSC,fname);
            else
                print(hAvPeriodHeatVsVsSC,fname,printDriver,dpi);
            end
        end
        if ~isempty(format)
            fname = [theme,'/VsFRVsSTDheat-',theme,'_',num2str(k-maxThetai,'%02d'),'.',format];
            if strcmp(format,'fig')
                saveas(hAvPeriodHeatVsVsSTD,fname);
            else
                print(hAvPeriodHeatVsVsSTD,fname,printDriver,dpi);
            end
        end
    end
end
function [binRate, binRateStd, VsID] = binVsFr(Vs,sVs,rate,VsEdge,sVsEdge, VsID)
    nbinVs = length(VsEdge)-1;
    nbinsVs = length(sVsEdge)-1;
    Vs = Vs(:);
    sVs = sVs(:);
    rate = rate(:);
    if nargin < 6
        [~,~,VsID] = histcounts(Vs,VsEdge);
    end
    [~,~,sVsID] = histcounts(Vs,sVsEdge);
    binRate = zeros(nbinVs,nbinsVs);
    binRateStd = zeros(nbinVs,nbinsVs);
    for i=1:nbinVs
        parfor j=1:nbinsVs
            pick = (VsID==i & sVsID == j);
            if sum(pick) > 0
                rt = rate(pick);
                binRate(i,j) = mean(rt);
                binRateStd(i,j) = std(rt);
            end
        end
    end
end
function heatVsStdFr(VsEdge,sVsEdge,rate,rateStd,ax1,ax2,sVslabel)
    nbinVs = length(VsEdge)-1;
    nbinsVs = length(sVsEdge)-1;
    max(rate(:))
    min(rate(:))
    axes(ax1);
    imagesc([1,nbinVs],[nbinsVs,1],rate');
    colormap('hot');
    colorbar;

    dVsEdge = VsEdge(2) - VsEdge(1);
    xtick = [1,ceil(nbinVs/2),nbinVs];
    xticklabel = num2str(VsEdge(1) + (xtick'-1)*dVsEdge);
    xtick = xtick-0.5;
    set(ax1,'XTick',xtick,'XTickLabel',xticklabel);
    xlabel(ax1,'Vs');
    title('FR');

    dsVsEdge = sVsEdge(2) - sVsEdge(1);
    ytick = [1,ceil(nbinVs/2),nbinVs];
    yticklabel = num2str(sVsEdge(1) + (ytick'-1)*dsVsEdge);
    ytick = ytick-0.5;
    set(ax1,'YTick',ytick,'YTickLabel',flipud(yticklabel));
    ylabel(ax1,sVslabel);

    axes(ax2);
    imagesc([1,nbinVs],[nbinsVs,1],rateStd');
    colormap('hot');
    colorbar;

    set(ax2,'XTick',xtick,'XTickLabel',xticklabel);
    xlabel(ax2,'Vs');
    set(ax2,'YTick',ytick,'YTickLabel',flipud(yticklabel));
    ylabel(ax2,sVslabel);
    title('\sigma(FR)');
end
