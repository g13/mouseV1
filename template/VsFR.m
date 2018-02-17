function VsFR(lgnfile,theme,ntheta,contrastLevel,ndperiod,loadData,format,dir)
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
    
    %set(0,'defaultAxesFontSize',6);
    %set(0,'defaultTextFontSize',8);    
    load(lgnfile,'p');
    epick = 1:p.nv1e;
    ipick = (p.nv1e+1):p.nv1;
    ntotal = p.nv1;
    dir = [dir,'/'];

    AvPeriodTransfer = true;
    NearInstantTransfer = false;
    explore = false;
    heat = true;
    transfer = true;
    meanOnly = true;

    % edges for Average ove period
    AvPVsEdgeE = {[0.05:0.01:0.35],[0.15:0.01:0.45],[0.35:0.025:0.85],[0.5:0.025:1.1]};
    AvPVsStdEdgeE = {[0.25:0.01:0.4],[0.3:0.01:0.46],[0.25:0.025:0.55],[0.25:0.025:0.575]};
    AvPVsScEdgeE = {[0:0.05:1.4],[0:0.05:1.2],[0:0.05:1.4],[0:0.05:1.0]};
    AvPVsEdgeI = {[0.15:0.01:0.55],[0.34:0.01:0.76],[0.7:0.025:1.4],[1.0:0.01:1.5]};
    AvPVsStdEdgeI ={[0.35:0.01:0.525],[0.4:0.01:0.54],[0.25:0.025:0.5],[0.225:0.01:0.425]};
    AvPVsScEdgeI = {[0:0.05:1.3],[0:0.05:0.6],[0:0.05:0.8],[0:0.01:0.4]};
    % edges for near instant
    NiVsEdgeE = 0:0.05:1.2;
    NiVsStdEdgeE = 0:0.05:0.6;
    NiVsScEdgeE = 0:0.1:2;
    NiVsEdgeI = 0:0.05:1.5;
    NiVsStdEdgeI = 0:0.05:0.6;
    NiVsScEdgeI = 0:0.1:2;
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
    tic;
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
                aC(i).spikes(j,:,:)  = aC(i).spikes(j,:,half);
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
                    aC(i).spikes(j,:,k)  = aC(i).spikes(j,iperiod,k);
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
    toc;
    % HSV color
    s0 = 0.5;
    h = (0:(ntheta-1))'./ntheta;
    s = s0 + (1-s0)*(0:(contrastLevel-1))./contrastLevel;
    v = 1.0;
    if AvPeriodTransfer
        if explore % for edges to histogram
            for k=1:ntheta
                hVs = figure;
                hFr = figure;
                for i=1:contrastLevel
                    if length(AvPVsEdgeE) > 1
                        VsEdgeE = AvPVsEdgeE{i};
                    else
                        VsEdgeE = AvPVsEdgeE;
                    end
                    if length(AvPVsStdEdgeE) > 1
                        VsStdEdgeE = AvPVsStdEdgeE{i};
                    else
                        VsStdEdgeE = AvPVsStdEdgeE;
                    end
                    if length(AvPVsScEdgeE) > 1
                        VsScEdgeE = AvPVsScEdgeE{i};
                    else
                        VsScEdgeE = AvPVsScEdgeE;
                    end
                    figure(hVs);
                    subplot(4,6,(i-1)*6+1)
                    histogram(Vs(epick,k,i));
                    if i==1, title('Exc Vs'), end;
                    ylabel(['C',num2str(i)]);
                    subplot(4,6,(i-1)*6+2)
                    histogram(VsSTD(epick,k,i));
                    if i==1, title('VsSTD'), end;
                    subplot(4,6,(i-1)*6+3)
                    histogram(VsSC(epick,k,i));
                    if i==1, title('VsSC'), end;

                    if length(AvPVsEdgeI) > 1
                        VsEdgeI = AvPVsEdgeI{i};
                    else
                        VsEdgeI = AvPVsEdgeI;
                    end
                    if length(AvPVsStdEdgeI) > 1
                        VsStdEdgeI = AvPVsStdEdgeI{i};
                    else
                        VsStdEdgeI = AvPVsStdEdgeI;
                    end
                    if length(AvPVsScEdgeI) > 1
                        VsScEdgeI = AvPVsScEdgeI{i};
                    else
                        VsScEdgeI = AvPVsScEdgeI;
                    end
                    subplot(4,6,(i-1)*6+4)
                    histogram(Vs(ipick,k,i));
                    if i==1, title('Inh Vs'), end;
                    subplot(4,6,(i-1)*6+5)
                    histogram(VsSTD(ipick,k,i));
                    if i==1, title('VsSTD'), end;
                    subplot(4,6,(i-1)*6+6)
                    histogram(VsSC(ipick,k,i));
                    if i==1, title('VsSC'), end;
                    figure(hFr);
                    subplot(2,4,i)
                    histogram(rate(epick,k,i));
                    if i==1, ylabel('Exc'),end;
                    title(['C',num2str(i)]);
                    subplot(2,4,4+i)
                    histogram(rate(ipick,k,i));
                    if i==1, ylabel('Inh'),end;
                end
                if ~isempty(format)
                    fname = [theme,'/AvPVsExplore-',theme,'_',num2str(k,'%02d'),'.',format];
                    if strcmp(format,'fig')
                        saveas(hVs,fname);
                    else
                        print(hVs,fname,printDriver,dpi);
                    end
                end
                if ~isempty(format)
                    fname = [theme,'/AvPFrExplore-',theme,'_',num2str(k,'%02d'),'.',format];
                    if strcmp(format,'fig')
                        saveas(hFr,fname);
                    else
                        print(hFr,fname,printDriver,dpi);
                    end
                end
            end
        end
        VsID_E = zeros(p.nv1e,contrastLevel,ntheta);
        VsID_I = zeros(p.nv1i,contrastLevel,ntheta);
        VsCID_E = zeros(p.nv1e,contrastLevel,ntheta);
        VsCID_I = zeros(p.nv1i,contrastLevel,ntheta);
        thetaRange = [1,7];
        contrastRange = [1,2,3,4];
        nc = length(contrastRange);
        if heat
        % Vs to Firing Rate transfer function averaged over period
            for ik=1:length(thetaRange)
                k = thetaRange(ik);
                hAvPHeatSTD = figure;
                hAvPHeatSC = figure;
                hAvPHeatStdSC = figure;
                for ii=1:nc
                    i = contrastRange(ii);
                    if length(AvPVsEdgeE) > 1
                        VsEdgeE = AvPVsEdgeE{i};
                    else
                        VsEdgeE = AvPVsEdgeE;
                    end
                    if length(AvPVsStdEdgeE) > 1
                        VsStdEdgeE = AvPVsStdEdgeE{i};
                    else
                        VsStdEdgeE = AvPVsStdEdgeE;
                    end
                    if length(AvPVsScEdgeE) > 1
                        VsScEdgeE = AvPVsScEdgeE{i};
                    else
                        VsScEdgeE = AvPVsScEdgeE;
                    end
                    figure(hAvPHeatSTD);
                    ax1 = subplot(nc,2,(ii-1)*2+1);
                    ax2 = [];
                        [binRateE,~,VsID] = binVsFr(Vs(epick,k,i),VsSTD(epick,k,i),rate(epick,k,i),VsEdgeE,VsStdEdgeE,[],false);
                        heatVsStdFr(VsEdgeE,VsStdEdgeE,binRateE,[],ax1,ax2);
                        if ii==1, title(ax1,'Exc FR'), end;
                        ylabel(ax1,['C',num2str(i),', \sigma(Vs)']);
                        if ii==nc, xlabel(ax1,'Vs'), end;
                        axis(ax1,'equal');
                    VsID_E(:,i,k) = VsID;

                    figure(hAvPHeatSC);
                    ax3 = subplot(nc,2,(ii-1)*2+1);
                    ax4 = [];
                        [binRateE,~,~] = binVsFr(Vs(epick,k,i),VsSC(epick,k,i),rate(epick,k,i),VsEdgeE,VsScEdgeE,VsID,false);
                        heatVsStdFr(VsEdgeE,VsScEdgeE,binRateE,[],ax3,ax4);
                        if ii==1, title(ax3,'Exc FR'), end;
                        ylabel(ax3,'Vs_{F1}');
                        if ii==nc, xlabel(ax3,'Vs'), end;

                    figure(hAvPHeatStdSC);
                    ax5 = subplot(nc,2,(ii-1)*2+1);
                    ax6 = [];
                        [binRateE,~,~] = binVsFr(VsSC(epick,k,i),VsSTD(epick,k,i),rate(epick,k,i),VsScEdgeE,VsStdEdgeE,[],false);
                        heatVsStdFr(VsScEdgeE,VsStdEdgeE,binRateE,[],ax5,ax6);
                        if ii==1, title(ax5,'Exc FR'), end;
                        ylabel(ax5,'\sigma(Vs)');
                        if ii==nc, xlabel(ax5,'Vs_{F1}'), end;
                    

                    if length(AvPVsEdgeI) > 1
                        VsEdgeI = AvPVsEdgeI{i};
                    else
                        VsEdgeI = AvPVsEdgeI;
                    end
                    if length(AvPVsStdEdgeI) > 1
                        VsStdEdgeI = AvPVsStdEdgeI{i};
                    else
                        VsStdEdgeI = AvPVsStdEdgeI;
                    end
                    if length(AvPVsScEdgeI) > 1
                        VsScEdgeI = AvPVsScEdgeI{i};
                    else
                        VsScEdgeI = AvPVsScEdgeI;
                    end
                    figure(hAvPHeatSTD);
                    ax1 = subplot(nc,2,(ii-1)*2+2);
                    ax2 = [];
                        [binRateI,~,VsID] = binVsFr(Vs(ipick,k,i),VsSTD(ipick,k,i),rate(ipick,k,i),VsEdgeI,VsStdEdgeI,[],false);
                        heatVsStdFr(VsEdgeI,VsStdEdgeI,binRateI,[],ax1,ax2);
                        if ii==1, title(ax1,'Inh FR'), end;
                        if ii==nc, xlabel(ax1,'Vs'), end;
                        axis(ax1,'equal');
                    VsID_I(:,i,k) = VsID;

                    figure(hAvPHeatSC);
                    ax3 = subplot(nc,2,(ii-1)*2+2);
                    ax4 = [];
                        [binRateI,~,~] = binVsFr(Vs(ipick,k,i),VsSC(ipick,k,i),rate(ipick,k,i),VsEdgeI,VsScEdgeI,VsID,false);
                        heatVsStdFr(VsEdgeI,VsScEdgeI,binRateI,[],ax3,ax4);
                        if ii==1, title(ax3,'Inh FR'), end;
                        if ii==nc, xlabel(ax3,'Vs'), end;

                    figure(hAvPHeatStdSC);
                    ax5 = subplot(nc,2,(ii-1)*2+2);
                    ax6 = [];
                        [binRateI,~,~] = binVsFr(VsSC(ipick,k,i),VsSTD(ipick,k,i),rate(ipick,k,i),VsScEdgeI,VsStdEdgeI,[],false);
                        heatVsStdFr(VsScEdgeI,VsStdEdgeI,binRateI,[],ax5,ax6);
                        if ii==1, title(ax5,'Inh FR'), end;
                        if ii==nc, xlabel(ax5,'Vs_{F1}'), end;
                end
                if ~isempty(format)
                    fname = [theme,'/AvPStdHeat-',theme,'_',num2str(k,'%02d')];
                    saveas(hAvPHeatSTD,fname);
                    print(hAvPHeatSTD,[fname,'.',format],printDriver,dpi);
                end
                if ~isempty(format)
                    fname = [theme,'/AvPScHeat-',theme,'_',num2str(k,'%02d')];
                    saveas(hAvPHeatSC,fname);
                    print(hAvPHeatSC,[fname,'.',format],printDriver,dpi);
                end
                if ~isempty(format)
                    fname = [theme,'/AvPStdScHeat-',theme,'_',num2str(k,'%02d')];
                    saveas(hAvPHeatStdSC,fname);
                    print(hAvPHeatStdSC,[fname,'.',format],printDriver,dpi);
                end
            end
        end

        if transfer
            hAvPVsFrTransferMeanOnly = figure;
            thetaRange = [1,7];
            for ik=1:length(thetaRange)
                k = thetaRange(ik);
                if ~meanOnly
                    hAvPVsFrTransfer = figure;
                end
                for i=1:contrastLevel
                    if length(AvPVsEdgeE) > 1
                        VsEdgeE = AvPVsEdgeE{i};
                    else
                        VsEdgeE = AvPVsEdgeE;
                    end
                    if length(AvPVsEdgeI) > 1
                        VsEdgeI = AvPVsEdgeI{i};
                    else
                        VsEdgeI = AvPVsEdgeI;
                    end
                    color = hsv2rgb(h(k),s(i),v);
                    if ~heat || sum(abs(i-contrastRange)<1e-14) == 0
                        [~,~,VsID_E(:,i,k)] = histcounts(Vs(epick,k,i),VsEdgeE);
                        [~,~,VsID_I(:,i,k)] = histcounts(Vs(ipick,k,i),VsEdgeI);
                    end
                    nbinsE = length(VsEdgeE)-1;
                    nbinsI = length(VsEdgeI)-1;
                    [~,VsCEdgeE,VsCID_E(:,i,k)] = histcounts(Vs(epick,k,i).*(1+VsSC(epick,k,i)),nbinsE);
                    [~,VsCEdgeI,VsCID_I(:,i,k)] = histcounts(Vs(ipick,k,i).*(1+VsSC(ipick,k,i)),nbinsI);

                    figure(hAvPVsFrTransferMeanOnly);
                    ax = subplot(2,2,1);
                    transferVsFr(VsID_E(:,i,k),rate(epick,k,i),VsEdgeE,ax,color);
                    if i==1 && ik==1 
                        hold on
                        xlabel('Vs');
                        ylabel('ExcFR (Hz)');
                    end

                    ax = subplot(2,2,2);
                    %transferVsFrMeanOnly(VsID_I(:,i,k),rate(ipick,k,i),VsEdgeI,ax,color);
                    transferVsFr(VsID_I(:,i,k),rate(ipick,k,i),VsEdgeI,ax,color);
                    if i==1 && ik==1 
                        hold on
                        xlabel('Vs');
                        ylabel('InhFR (Hz)');
                    end

                    ax = subplot(2,2,3);
                    transferVsFr(VsCID_E(:,i,k),rate(epick,k,i),VsCEdgeE,ax,color);
                    if i==1 && ik==1 
                        hold on
                        xlabel('V_{s}^{F0+F1}');
                        ylabel('ExcFR (Hz)');
                    end

                    ax = subplot(2,2,4);
                    transferVsFr(VsCID_I(:,i,k),rate(ipick,k,i),VsCEdgeI,ax,color);
                    if i==1 && ik==1 
                        hold on
                        xlabel('V_{s}^{F0+F1}');
                        ylabel('InhFR (Hz)');
                    end

                    if ~meanOnly
                        figure(hAvPVsFrTransfer);
                        ax = subplot(2,contrastLevel,i);
                            transferVsFr(VsID_E(:,i,k),rate(epick,k,i),VsEdgeE,ax);
                            title(ax,['C',num2str(i)]);
                            if i==1, ylabel(ax,'Exc FR'), end;
                        ax = subplot(2,contrastLevel,i+4);
                            transferVsFr(VsID_I(:,i,k),rate(ipick,k,i),VsEdgeI,ax);
                            if i==1, ylabel(ax,'Inh FR'), end;
                            if i==1, xlabel(ax,'Vs'), end;
                    end
                end
                if ~isempty(format) && ~meanOnly
                    fname = [theme,'/AvPTransferVsFR-',theme,'_',num2str(k,'%02d')];
                    saveas(hAvPVsFrTransfer,fname);
                    print(hAvPVsFrTransfer,[fname,'.',format],printDriver,dpi);
                end
            end
            if ~isempty(format)
                fname = [theme,'/AvPFRVsMeanOnly-',theme];
                saveas(hAvPVsFrTransferMeanOnly,fname);
                print(hAvPVsFrTransferMeanOnly,[fname,'.',format],printDriver,dpi);
            end
        end
    end
    periodSample = [1,maxPeriodi,2*maxPeriodi-1,3*maxPeriodi-2];
    v = zeros(ndperiod,1);
    v(periodSample) = [0.7,1.0,0.7,0.4];
    nPeriodSample = length(periodSample);
    if NearInstantTransfer
        VsEdgeE = NiVsEdgeE;
        VsStdEdgeE =NiVsStdEdgeE;
        VsScEdgeE = NiVsScEdgeE;
        VsEdgeI = NiVsEdgeI;
        VsStdEdgeI = NiVsStdEdgeI;
        VsScEdgeI = NiVsScEdgeI;
        if explore % for edges to histogram
            for k=1:ntheta
                for j=periodSample
                    hVs = figure;
                    hFr = figure;
                    for i=1:contrastLevel
                        figure(hVs);
                        subplot(contrastLevel,4,(i-1)*4+1)
                        histogram(aC(i).Vs(epick,j,k));
                        if i==1, title(['Exc Vs ',num2str(j),'/',num2str(ndperiod),'period']), end;
                        ylabel(['C',num2str(i)]);

                        subplot(contrastLevel,4,(i-1)*4+2)
                        histogram(aC(i).VsSTD(epick,j,k));
                        if i==1, title('VsSTD'), end;

                        subplot(contrastLevel,4,(i-1)*4+3)
                        histogram(aC(i).Vs(ipick,j,k));
                        if i==1, title('Inh Vs'), end;

                        subplot(contrastLevel,4,(i-1)*4+4)
                        histogram(aC(i).VsSTD(ipick,j,k));
                        if i==1, title('VsSTD'), end;

                        figure(hFr);
                        subplot(2,contrastLevel,i)
                        histogram(aC(i).spikes(epick,j,k));
                        if i==1, ylabel(['Exc FR ',num2str(j),'/',num2str(ndperiod),'period']),end;
                        title(['C',num2str(i)]);
                        subplot(2,contrastLevel,contrastLevel+i)
                        histogram(aC(i).spikes(ipick,j,k));
                        if i==1, ylabel('Inh'),end;
                    end
                    if ~isempty(format)
                        fname = [theme,'/NiVsExplore-',theme,'_',num2str(k,'%02d'),'-',num2str(j)];
                        saveas(hVs,fname);
                        print(hVs,[fname,'.',format],printDriver,dpi);
                    end
                    if ~isempty(format)
                        fname = [theme,'/NiFrExplore-',theme,'_',num2str(k,'%02d'),'-',num2str(j)];
                        saveas(hFr,fname);
                        print(hFr,[fname,'.',format],printDriver,dpi);
                    end
                end
            end
        end
        VsID_E = zeros(p.nv1e,contrastLevel,nPeriodSample,ntheta);
        VsID_I = zeros(p.nv1i,contrastLevel,nPeriodSample,ntheta);
        if heat
        % Vs to Firing Rate transfer function averaged over period
            for k=1:ntheta
                for jj=1:nPeriodSample
                    j = periodSample(jj);
                    hNiHeatVsVsSTD = figure;
                    for i=1:contrastLevel
                        figure(hNiHeatVsVsSTD);
                        ax1 = subplot(4,4,(i-1)*4+1);
                        ax2 = subplot(4,4,(i-1)*4+2);
                            [binRateE,binRateStdE,VsID] = binVsFr(aC(i).Vs(epick,j,k),aC(i).VsSTD(epick,j,k),aC(i).spikes(epick,j,k),VsEdgeE,VsStdEdgeE);
                            heatVsStdFr(VsEdgeE,VsStdEdgeE,binRateE,binRateStdE,ax1,ax2);
                            if i==1, title(ax1,['Exc FR',num2str(j),'/',num2str(ndperiod),'period']), end;
                            if i==1, title(ax2,'\sigma(FR)'), end;
                            ylabel(ax1,['C',num2str(i),', \sigma(Vs)']);
                            if i==contrastLevel, xlabel(ax1,'Vs'), end;
                        VsID_E(:,i,jj,k) = VsID;

                        figure(hNiHeatVsVsSTD);
                        ax1 = subplot(4,4,(i-1)*4+3);
                        ax2 = subplot(4,4,(i-1)*4+4);
                            [binRateI,binRateStdI,VsID] = binVsFr(aC(i).Vs(ipick,j,k),aC(i).VsSTD(ipick,j,k),aC(i).spikes(ipick,j,k),VsEdgeI,VsStdEdgeI);
                            heatVsStdFr(VsEdgeI,VsStdEdgeI,binRateI,binRateStdI,ax1,ax2);
                            if i==1, title(ax1,'Inh FR'), end;
                            if i==1, title(ax2,'\sigma(FR)'), end;
                        VsID_I(:,i,jj,k) = VsID;

                    end
                    if ~isempty(format)
                        fname = [theme,'/NiVsFRVsSTDheat-',theme,'_',num2str(k,'%02d'),'-',num2str(j)];
                        saveas(hNiHeatVsVsSTD,fname);
                        print(hNiHeatVsVsSTD,[fname,'.',format],printDriver,dpi);
                    end
                end
            end
        end

        if transfer
            hNiVsFrTransferMeanOnly = figure;
            for k=1:ntheta
                for jj=1:nPeriodSample
                    j = periodSample(jj);
                    if ~meanOnly
                        hNiVsFrTransfer = figure;
                    end
                    for i=1:contrastLevel
                        if ~heat
                            [~,~,VsID_E(:,i,jj,k)] = histcounts(aC(i).Vs(epick,j,k),VsEdgeE);
                            [~,~,VsID_I(:,i,jj,k)] = histcounts(aC(i).Vs(ipick,j,k),VsEdgeI);
                        end
                        figure(hNiVsFrTransferMeanOnly);
                        color = hsv2rgb(h(k),s(i),v(j));
                        ax = subplot(1,2,1);
                        transferVsFrMeanOnly(VsID_E(:,i,jj,k),aC(i).spikes(epick,j,k),VsEdgeE,ax,color);
                        if i==1 && k==1 && jj==1
                            hold on
                            xlabel('Vs');
                            ylabel('ExcFR (Hz)');
                        end

                        ax = subplot(1,2,2);
                        transferVsFrMeanOnly(VsID_I(:,i,jj,k),aC(i).spikes(ipick,j,k),VsEdgeI,ax,color);
                        if i==1 && k==1 && jj==1
                            hold on
                            xlabel('Vs');
                            ylabel('InhFR (Hz)');
                        end
                        if ~meanOnly
                            figure(hNiVsFrTransfer);
                            ax = subplot(2,4,i);
                                transferVsFr(VsID_E(:,i,jj,k),aC(i).spikes(epick,j,k),VsEdgeE,ax);
                                title(ax,['C',num2str(i)]);
                                if i==1, ylabel(ax,['Exc FR ', num2str(j),'/',num2str(ndperiod),'period']), end;
                            ax = subplot(2,4,i+4);
                                transferVsFr(VsID_I(:,i,jj,k),aC(i).spikes(ipick,j,k),VsEdgeI,ax);
                                if i==1, ylabel(ax,'Inh FR'), end;
                                if i==1, xlabel(ax,'Vs'), end;
                        end
                    end
                    if ~isempty(format) && ~meanOnly
                        fname = [theme,'/NiTransferVsFR-',theme,'_',num2str(k,'%02d'),'-',num2str(j),'.',format];
                        saveas(hNiVsFrTransfer,fname);
                        print(hNiVsFrTransfer,[fname,'.',format],printDriver,dpi);
                    end
                end
            end
            if ~isempty(format)
                fname = [theme,'/NiFRVsMeanOnly-',theme];
                saveas(hNiVsFrTransferMeanOnly,fname);
                print(hNiVsFrTransferMeanOnly,[fname,'.',format],printDriver,dpi);
            end
        end
    end
end
function [binRate, binRateStd, VsID] = binVsFr(Vs,sVs,rate,VsEdge,sVsEdge, VsID, binStd)
    nbinVs = length(VsEdge)-1;
    nbinsVs = length(sVsEdge)-1;
    Vs = Vs(:);
    sVs = sVs(:);
    rate = rate(:);
    if nargin < 7
        binStd = true;
    end
    if nargin < 6 || isempty(VsID)
        [~,~,VsID] = histcounts(Vs,VsEdge);
    end
    [~,~,sVsID] = histcounts(sVs,sVsEdge);
    binRate = zeros(nbinVs,nbinsVs);
    if binStd
        binRateStd = zeros(nbinVs,nbinsVs);
    else
        binRateStd = [];
    end
    for i=1:nbinVs
        for j=1:nbinsVs
            pick = (VsID==i & sVsID == j);
            if sum(pick) > 0
                rt = rate(pick);
                binRate(i,j) = mean(rt);
                if binStd
                    binRateStd(i,j) = std(rt);
                end
            end
        end
    end
end
function heatVsStdFr(VsEdge,sVsEdge,rate,rateStd,ax1,ax2)
    nbinVs = length(VsEdge)-1;
    nbinsVs = length(sVsEdge)-1;
    imagesc([1,nbinVs],[nbinsVs,1],rate','Parent',ax1);
    colormap('hot');
    colorbar('peer',ax1);
    box off

    dVsEdge = VsEdge(2) - VsEdge(1);
    %xtick = [1,ceil(nbinVs/2),nbinVs];
    xtick = get(ax1,'xtick');
    xticklabel = num2str(VsEdge(1) + (xtick'-1)*dVsEdge);
    xtick = xtick-0.5;
    set(ax1,'XTick',xtick,'XTickLabel',xticklabel);

    dsVsEdge = sVsEdge(2) - sVsEdge(1);
    %ytick = [1,ceil(nbinsVs/2),nbinsVs];
    ytick = get(ax1,'ytick');
    yticklabel = num2str(sVsEdge(1) + (ytick'-1)*dsVsEdge);
    ytick = ytick-0.5;
    set(ax1,'YTick',ytick,'YTickLabel',flipud(yticklabel));
    xlim([0,nbinVs+1])
    ylim([0,nbinsVs+1])

    if nargin ==6 && ~isempty(ax2)
        imagesc([1,nbinVs],[nbinsVs,1],rateStd','Parent',ax2);
        colormap('hot');
        colorbar('peer',ax2);

        set(ax2,'XTick',xtick,'XTickLabel',xticklabel);
        set(ax2,'YTick',ytick,'YTickLabel',flipud(yticklabel));
    end
end
function transferVsFr(VsID,rate,VsEdge,ax,color);
    if nargin<5 
        color = [];
    end
    dVsEdge = VsEdge(2) - VsEdge(1);
    nbins = length(VsEdge)-1;
    mRate = zeros(nbins,1);
    uRate = zeros(nbins,1);
    lRate = zeros(nbins,1);
    for i = 1:nbins
        mRate(i) = mean(rate(VsID==i)); 
        lRate(i) = mRate(i)-quantile(rate(VsID==i),0.25); 
        uRate(i) = quantile(rate(VsID==i),0.75)-mRate(i); 
    end
    if isempty(color)
        errorbar(ax,VsEdge(1:nbins)+dVsEdge,mRate,lRate,uRate);
    else
        errorbar(ax,VsEdge(1:nbins)+dVsEdge,mRate,lRate,uRate,'Color',color);
    end
end
function transferVsFrMeanOnly(VsID,rate,VsEdge,ax,color);
    dVsEdge = VsEdge(2) - VsEdge(1);
    nbins = length(VsEdge)-1;
    mRate = zeros(nbins,1);
    for i = 1:nbins
        mRate(i) = mean(rate(VsID==i)); 
    end
    plot(ax,VsEdge(1:nbins)+dVsEdge,mRate,'LineWidth',1.5,'Color',color);
end
