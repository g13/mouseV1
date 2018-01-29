function simAssisVeff(theme,loadData,fitReady,format,dir,npool)
    if nargin <6
        npool = 6;
        if nargin < 5
            dir = '.';
            if nargin < 4
                format = 'fig';
                if nargin < 3
                    fitReady = false;
                    if nargin < 2
                        loadData = true;
                    end
                end
            end
        end
    end
    dir = [dir,'/'];
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r300';
    end
    Vs_0 = 0:0.05:0.5;
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'DefaultErrorbarLineWidth', 2);
    set(0, 'DefaultAxesFontSize', 16);
    n = 10800;
    ntheta = 12;
    dtheta = 180/ntheta;
    ne = 8460;
    contrastLevel = 4;
    contrast = 1:contrastLevel;
    ndperiod = 25;
    vE = 2.8;
    vI = -0.4;
    vL = 0;
    EgL = 50;
    IgL = 70;
    conLabel = cell(contrastLevel,1);
    if ~loadData
        for i=1:contrastLevel
            disp(['loading c',num2str(i)]);
            eps = 125*2^(i-1);
            conLabel(i) = {[num2str(eps/10,'%.1f'),'%']};
            DIR = [theme,'/',num2str(eps,'%04d')];
            [nP(i),tC(i)] = readSimple(DIR,ntheta,n,ndperiod);
            % tC [n,2*ntheta+1] //frate,stifr [n, 2*ntheta]
            % nP [n,1]
        end
        save([theme,'Simple.mat'],'tC','nP');
    else 
        load([theme,'Simple.mat'],'tC','nP');
    end
    if ~fitReady
        Veff = zeros(n,2*ntheta,contrastLevel);
        prA = zeros(n,contrastLevel);
        for i=1:contrastLevel
            Veff(:,:,i) = tC(i).Veff(:,1:2*ntheta).*tC(i).Veffsc(:,1:2*ntheta);
            prA(:,i) = nP(i).priA;
        end
        [Veffmax, prefAngle, sigfsq2, lifted] = fitVeffTC(Veff,prA,contrastLevel,n,ntheta,npool);
        save([theme,'VeffFit.mat'],'Veffmax','prefAngle','sigfsq2','lifted');
    else
        load([theme,'VeffFit.mat']);
    end
    nVs_0 = length(Vs_0);
    figure
    %  Veff E        % Veff I     
    %  norm Veff E   % norm Veff I
    
    evenOdd = mod(ntheta,2);
    relTheta = (-(ntheta+evenOdd)/2:(ntheta+evenOdd)/2)*dtheta;
    TCVeff_f0 = zeros(ntheta+1,n,contrastLevel)-1;
    TCVeff_f0_f1 = TCVeff_f0;
    %TCrate = TCVeff_f0;
    LineScheme = {':','-.','--','-'};
    epick = 1:ne;
    ipick = (ne+1):n;
    for i = 1:contrastLevel
        for k=1:n    
            ipA = round(nP(i).priA(k)*180/pi/dtheta)+1;
            if ipA > ntheta
                ipA = ipA-ntheta;
            end
            if ipA < ntheta/2+1
                half = [(ipA+ntheta/2):ntheta,1:ipA+ntheta/2];
            else 
                half = [ipA-ntheta/2:ntheta,1:(ipA-ntheta/2)];
            end
            %rho = [tC(i).frate(k,:),tC(i).frate(k,1)];
            Veff_f0 = [tC(i).Veff(k,:),tC(i).Veff(k,1)];
            Veff_f0_f1 = Veff_f0.*(1+[tC(i).Veffsc(k,:),tC(i).Veffsc(k,1)]);
            %TCrate(:,k,i) = rho(half);
            TCVeff_f0(:,k,i) = Veff_f0(half);
            TCVeff_f0_f1(:,k,i) = Veff_f0_f1(half);
        end
        subplot(1,2,1)
            hold on
            plot(relTheta,mean(TCVeff_f0(:,epick,i),2),'r','LineStyle',LineScheme{i});
            plot(relTheta,mean(TCVeff_f0(:,ipick,i),2),'b','LineStyle',LineScheme{i});
        subplot(1,2,2)
            hold on 
            plot(relTheta,mean(TCVeff_f0_f1(:,epick,i),2),'r','LineStyle',LineScheme{i});
            plot(relTheta,mean(TCVeff_f0_f1(:,ipick,i),2),'b','LineStyle',LineScheme{i});
    end
    subplot(1,2,1)
    legend({'exc','inh'});
    title('Veff F0');
    xlabel('\theta');
    ylim([0,inf]);
    subplot(1,2,2)
    title('Veff F0+F1');
    xlabel('\theta');
    ylim([0,inf]);
    filename = [dir,'VeffTC'];
    
    if strcmp(format,'fig')
        saveas(gcf,[filename,'.',format]);
    else
        print(gcf,[filename,'.',format],printDriver,dpi);
    end
    figure
    %  OSI Vs F0+F1(Vs_0)       %  OSI Vs F0(Vs_0)      % sudo f0_f1(Vs_0) % sudo f0(Vs_0)
    %  OSI Vs F0+F1(constrast)  %  OSI Vs F0(contrast)  % sudo mean(contrast)
    Vs_f0_f1= zeros(n,2,contrastLevel);
    Vs_f0 = zeros(n,2,contrastLevel);
    Vs_f0_f1_s = zeros(n,2,contrastLevel);
    Vs_f0_s = zeros(n,2,contrastLevel);
    Vs_0_mat = ones(n,1) * Vs_0;
    for i=1:contrastLevel
        for j=1:n
            itheta = nP(i).indpo(j);
            Vs_f0_f1(j,1,i) = tC(i).Veff(j,itheta)*(1+tC(i).Veffsc(j,itheta));
            Vs_f0(j,1,i) = tC(i).Veff(j,itheta);
            if j<=ne
                Vs_f0_f1_s(j,1,i) = ((tC(i).gE(j,itheta)*(1+tC(i).gEsc(j,itheta))+tC(i).gLGN(j,itheta)*(1+tC(i).gLGNsc(j,itheta)))*vE + tC(i).gI(j,itheta)*vI + EgL*vL)/tC(i).gtot(j,itheta);
                Vs_f0_s(j,1,i) = ((tC(i).gE(j,itheta)+tC(i).gLGN(j,itheta))*vE + tC(i).gI(j,itheta)*vI + EgL*vL)/tC(i).gtot(j,itheta);
            else
                Vs_f0_f1_s(j,1,i) = ((tC(i).gE(j,itheta)*(1+tC(i).gEsc(j,itheta))+tC(i).gLGN(j,itheta)*(1+tC(i).gLGNsc(j,itheta)))*vE + tC(i).gI(j,itheta)*vI + IgL*vL)/tC(i).gtot(j,itheta);
                Vs_f0_s(j,1,i) = ((tC(i).gE(j,itheta)+tC(i).gLGN(j,itheta))*vE + tC(i).gI(j,itheta)*vI + IgL*vL)/tC(i).gtot(j,itheta);
            end
            itheta = nP(i).indoo(j);
            Vs_f0_f1(j,2,i) = tC(i).Veff(j,itheta)*(1+tC(i).Veffsc(j,itheta));
            Vs_f0(j,2,i) = tC(i).Veff(j,itheta);
            if j<=ne
                Vs_f0_f1_s(j,2,i) = ((tC(i).gE(j,itheta)*(1+tC(i).gEsc(j,itheta))+tC(i).gLGN(j,itheta)*(1+tC(i).gLGNsc(j,itheta)))*vE + tC(i).gI(j,itheta)*vI + EgL*vL)/tC(i).gtot(j,itheta);
                Vs_f0_s(j,2,i) = ((tC(i).gE(j,itheta)+tC(i).gLGN(j,itheta))*vE + tC(i).gI(j,itheta)*vI + EgL*vL)/tC(i).gtot(j,itheta);
            else
                Vs_f0_f1_s(j,2,i) = ((tC(i).gE(j,itheta)*(1+tC(i).gEsc(j,itheta))+tC(i).gLGN(j,itheta)*(1+tC(i).gLGNsc(j,itheta)))*vE + tC(i).gI(j,itheta)*vI + IgL*vL)/tC(i).gtot(j,itheta);
                Vs_f0_s(j,2,i) = ((tC(i).gE(j,itheta)+tC(i).gLGN(j,itheta))*vE + tC(i).gI(j,itheta)*vI + IgL*vL)/tC(i).gtot(j,itheta);
            end
        end
    end
    subplot(2,4,1)
        hold on
        osi_f0_f1_h = osi(Vs_f0_f1(:,1,4)*ones(1,nVs_0)-Vs_0_mat,Vs_f0_f1(:,2,4)*ones(1,nVs_0)-Vs_0_mat);
        osi_f0_f1_l = osi(Vs_f0_f1(:,1,2)*ones(1,nVs_0)-Vs_0_mat,Vs_f0_f1(:,2,2)*ones(1,nVs_0)-Vs_0_mat);
        %errorbar(Vs_0,mean(osi_f0_f1_h(1:ne,:)),std(osi_f0_f1_h(1:ne,:)),'-r');
        %errorbar(Vs_0,mean(osi_f0_f1_l(1:ne,:)),std(osi_f0_f1_l(1:ne,:)),':r');
        %errorbar(Vs_0,mean(osi_f0_f1_h((ne+1):n,:)),std(osi_f0_f1_h((ne+1):n,:)),'-b');
        %errorbar(Vs_0,mean(osi_f0_f1_l((ne+1):n,:)),std(osi_f0_f1_l((ne+1):n,:)),':b');
        plot(Vs_0,mean(osi_f0_f1_h(1:ne,:)),'-r');
        plot(Vs_0,mean(osi_f0_f1_l(1:ne,:)),':r');
        plot(Vs_0,mean(osi_f0_f1_h((ne+1):n,:)),'-b');
        plot(Vs_0,mean(osi_f0_f1_l((ne+1):n,:)),':b');
        legend({'high exc','low exc','high inh','low inh'},'FontSize',10);
        title('high as 100% (4th), low as 25% (2ed)','FontSize',10); 
        ylim([0,inf]);
        xlabel('Vs_0');
        ylabel('True OSI Vs F0+F1 - Vs_0');

    subplot(2,4,2)
        hold on
        osi_f0_h = osi(Vs_f0(:,1,4)*ones(1,nVs_0)-Vs_0_mat,Vs_f0(:,2,4)*ones(1,nVs_0)-Vs_0_mat);
        osi_f0_l = osi(Vs_f0(:,1,2)*ones(1,nVs_0)-Vs_0_mat,Vs_f0(:,2,2)*ones(1,nVs_0)-Vs_0_mat);
        %errorbar(Vs_0,mean(osi_f0_h(1:ne,:)),std(osi_f0_h(1:ne,:)),'-r');
        %errorbar(Vs_0,mean(osi_f0_l(1:ne,:)),std(osi_f0_l(1:ne,:)),':r');
        %errorbar(Vs_0,mean(osi_f0_h((ne+1):n,:)),std(osi_f0_h((ne+1):n,:)),'-b');
        %errorbar(Vs_0,mean(osi_f0_l((ne+1):n,:)),std(osi_f0_l((ne+1):n,:)),':b');
        plot(Vs_0,mean(osi_f0_h(1:ne,:)),'-r');
        plot(Vs_0,mean(osi_f0_l(1:ne,:)),':r');
        plot(Vs_0,mean(osi_f0_h((ne+1):n,:)),'-b');
        plot(Vs_0,mean(osi_f0_l((ne+1):n,:)),':b');
        ylim([0,inf]);
        xlabel('Vs_0');
        ylabel('True OSI Vs F0 - Vs_0');

    subplot(2,4,3)
        hold on
        osi_f0_f1_s_h = osi(Vs_f0_f1_s(:,1,4)*ones(1,nVs_0)-Vs_0_mat,Vs_f0_f1_s(:,2,4)*ones(1,nVs_0)-Vs_0_mat);
        osi_f0_f1_s_l = osi(Vs_f0_f1_s(:,1,2)*ones(1,nVs_0)-Vs_0_mat,Vs_f0_f1_s(:,2,2)*ones(1,nVs_0)-Vs_0_mat);
        %errorbar(Vs_0,mean(osi_f0_f1_s_h(1:ne,:)),std(osi_f0_f1_s_h(1:ne,:)),'-r');
        %errorbar(Vs_0,mean(osi_f0_f1_s_l(1:ne,:)),std(osi_f0_f1_s_l(1:ne,:)),':r');
        %errorbar(Vs_0,mean(osi_f0_f1_s_h((ne+1):n,:)),std(osi_f0_f1_s_h((ne+1):n,:)),'-b');
        %errorbar(Vs_0,mean(osi_f0_f1_s_l((ne+1):n,:)),std(osi_f0_f1_s_l((ne+1):n,:)),':b');
        plot(Vs_0,mean(osi_f0_f1_s_h(1:ne,:)),'-r');
        plot(Vs_0,mean(osi_f0_f1_s_l(1:ne,:)),':r');
        plot(Vs_0,mean(osi_f0_f1_s_h((ne+1):n,:)),'-b');
        plot(Vs_0,mean(osi_f0_f1_s_l((ne+1):n,:)),':b');
        ylim([0,inf]);
        xlabel('Vs_0');
        ylabel('Pseudo OSI Vs F0+F1 - Vs_0');

    subplot(2,4,4)
        hold on
        osi_f0_s_h = osi(Vs_f0_s(:,1,4)*ones(1,nVs_0)-Vs_0_mat,Vs_f0_s(:,2,4)*ones(1,nVs_0)-Vs_0_mat);
        osi_f0_s_l = osi(Vs_f0_s(:,1,2)*ones(1,nVs_0)-Vs_0_mat,Vs_f0_s(:,2,2)*ones(1,nVs_0)-Vs_0_mat);
        %errorbar(Vs_0,mean(osi_f0_s_h(1:ne,:)),std(osi_f0_s_h(1:ne,:)),'-r');
        %errorbar(Vs_0,mean(osi_f0_s_l(1:ne,:)),std(osi_f0_s_l(1:ne,:)),':r');
        %errorbar(Vs_0,mean(osi_f0_s_h((ne+1):n,:)),std(osi_f0_s_h((ne+1):n,:)),'-b');
        %errorbar(Vs_0,mean(osi_f0_s_l((ne+1):n,:)),std(osi_f0_s_l((ne+1):n,:)),':b');
        plot(Vs_0,mean(osi_f0_s_h(1:ne,:)),'-r');
        plot(Vs_0,mean(osi_f0_s_l(1:ne,:)),':r');
        plot(Vs_0,mean(osi_f0_s_h((ne+1):n,:)),'-b');
        plot(Vs_0,mean(osi_f0_s_l((ne+1):n,:)),':b');
        ylim([0,inf]);
        xlabel('Vs_0');
        ylabel('Pseudo OSI Vs F0 - Vs_0');

    subplot(2,4,5)
        hold on
        osi_f0_f1 = osi(squeeze(Vs_f0_f1(:,1,:)),squeeze(Vs_f0_f1(:,2,:)));
        errorbar(contrast,mean(osi_f0_f1(1:ne,:)),std(osi_f0_f1(1:ne,:)),':r');
        errorbar(contrast,mean(osi_f0_f1((ne+1):n,:)),std(osi_f0_f1((ne+1):n,:)),':b');
        xlim([0,5])
        set(gca,'XTick',[1,2,3,4]);
        set(gca,'XTickLabel',{'12.5','25','50','100'});
        xlabel('contrast %');
        ylabel('True OSI Vs F0+F1');

    subplot(2,4,6)
        hold on
        osi_f0 = osi(squeeze(Vs_f0(:,1,:)),squeeze(Vs_f0(:,2,:)));
        errorbar(contrast,mean(osi_f0(1:ne,:)),std(osi_f0(1:ne,:)),':r');
        errorbar(contrast,mean(osi_f0((ne+1):n,:)),std(osi_f0((ne+1):n,:)),':b');
        xlim([0,5])
        set(gca,'XTick',[1,2,3,4]);
        set(gca,'XTickLabel',{'12.5','25','50','100'});
        xlabel('contrast %');
        ylabel('True OSI Vs F0');

    subplot(2,4,7)
        hold on
        osi_f0_f1_s = osi(squeeze(Vs_f0_f1_s(:,1,:)),squeeze(Vs_f0_f1_s(:,2,:)));
        errorbar(contrast,mean(osi_f0_f1_s(1:ne,:)),std(osi_f0_f1_s(1:ne,:)),':r');
        errorbar(contrast,mean(osi_f0_f1_s((ne+1):n,:)),std(osi_f0_f1_s((ne+1):n,:)),':b');
        xlim([0,5])
        set(gca,'XTick',[1,2,3,4]);
        set(gca,'XTickLabel',{'12.5','25','50','100'});
        xlabel('contrast %');
        ylabel('Pseudo OSI Vs F0+F1');

    subplot(2,4,8)
        hold on
        osi_f0_s = osi(squeeze(Vs_f0_s(:,1,:)),squeeze(Vs_f0_s(:,2,:)));
        errorbar(contrast,mean(osi_f0_s(1:ne,:)),std(osi_f0_s(1:ne,:)),':r');
        errorbar(contrast,mean(osi_f0_s((ne+1):n,:)),std(osi_f0_s((ne+1):n,:)),':b');
        xlim([0,5])
        set(gca,'XTick',[1,2,3,4]);
        set(gca,'XTickLabel',{'12.5','25','50','100'});
        xlabel('contrast %');
        ylabel('Pseudo OSI Vs F0');
        filename = [dir,'OSIVeff'];
        if strcmp(format,'fig')
            saveas(gcf,[filename,'.',format]);
        else
            print(gcf,[filename,'.',format],printDriver,dpi);
        end

    figure;
    p2 = contrastLevel; p1 = contrastLevel-2;
    thres = 0.0;
    width = sigfsq2*sqrt(log(2));
    %pick = nP(contrastLevel).ei>0.5 & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(p1).pkrate>thres & nP(p2).pkrate>thres;
    pick = nP(contrastLevel).ei>0.5 & sum(tC(p1).frate>0,2)>=ntheta;
    subplot(2,4,1); hold on;
        edges = 0:5:90;
        pickedWidth = width(pick,p1);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','r');
        title('25%');
        xlabel('half width (\theta)')
        ylabel('% exc neurons')
    xlim([0,90]);
    subplot(2,4,2); hold on;
        edges = 0:5:90;
        pickedWidth = width(pick,p2);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','r');
        title('100%');
        xlabel('half width (\theta)')
        ylabel('% exc neurons')
    xlim([0,90]);

    hExcWidth = subplot(2,2,2);
        hold on
        idTick = 8;
        ctrs = cell(2,1);
        dTick = 30;
        tickLim = [0,90];
        if sum(pick) > 0
            pair1 = width(pick,p1);
            pair2 = width(pick,p2);
            [tick,n0] = autoAxis(tickLim(1),tickLim(2),round(tickLim(2)/dTick),tickLim);
            tickLabel = num2str(tick');
            lctrsx = (n0-1)*idTick+1;
            lctrsy = lctrsx;
            ctrs{1} = linspace(tick(1),tick(end),lctrsx);
            ctrs{2} = ctrs{1};
            tickPosY = 0.5:idTick:(lctrsy-1+0.5);
            tickPosX = 0.5:idTick:(lctrsx-1+0.5);
            
            dataPair = [pair1, pair2];
            denPair = hist3(dataPair,ctrs);
            maxDen = max(max(denPair));
            denPair = denPair/maxDen;
            denPair = denPair(1:(lctrsx-1),1:(lctrsy-1));
            imagesc([1,lctrsx-1],[1,lctrsy-1],denPair');
            plot(0.5:lctrsx+0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
            title([num2str(sum(pick)/sum(nP(contrastLevel).ei>0.5)*100,'%.1f'),'% qualified neurons']);
            xlabel('half width(25%)')
            ylabel('half width(100%)')
            daspect([1,1,1]);
            
            set(gca,'YTickLabel',tickLabel,'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
            colormap(hExcWidth,redOnly);
        end

    pick = nP(contrastLevel).ei<0.5 & sum(tC(p1).frate>0,2)>=ntheta;
    %pickedWidth = width(nP(p1).ei<0.5 & nP(p1).pkrate> nP(1).br & nP(p1).pkrate>thres,p1);
    subplot(2,4,5); hold on;
        edges = 0:5:90;
        pickedWidth = width(pick,p1);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','b');
        title('25%');
        xlabel('half width (\theta)')
        ylabel('% inh neurons')
        xlim([0,90]);
    subplot(2,4,6); hold on;
        edges = 0:5:90;
        pickedWidth = width(pick,p2);
        histogram(pickedWidth,edges,'Normalization','probability','FaceColor','b');
        title('100%');
        xlabel('half width (\theta)')
        ylabel('% inh neurons')
        xlim([0,90]);

    hInhWidth = subplot(2,2,4); 
        hold on
        idTick = 8;
        ctrs = cell(2,1);
        dTick = 30;
        tickLim = [0,90];
        if sum(pick) > 0
            pair1 = width(pick,p1);
            pair2 = width(pick,p2);
            [tick,n0] = autoAxis(tickLim(1),tickLim(2),round(tickLim(2)/dTick),tickLim);
            tickLabel = num2str(tick');
            lctrsx = (n0-1)*idTick+1;
            lctrsy = lctrsx;
            ctrs{1} = linspace(tick(1),tick(end),lctrsx);
            ctrs{2} = ctrs{1};
            tickPosY = 0.5:idTick:(lctrsy-1+0.5);
            tickPosX = 0.5:idTick:(lctrsx-1+0.5);
            
            dataPair = [pair1, pair2];
            denPair = hist3(dataPair,ctrs);
            maxDen = max(max(denPair));
            denPair = denPair/maxDen;
            denPair = denPair(1:(lctrsx-1),1:(lctrsy-1));
            imagesc([1,lctrsx-1],[1,lctrsy-1],denPair');
            plot(0.5:lctrsx+0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
            title([num2str(sum(pick)/sum(nP(contrastLevel).ei<0.5)*100,'%.1f'),'% qualified neurons']);
            xlabel('half width(25%)')
            ylabel('half width(100%)')
            daspect([1,1,1]);
            
            set(gca,'YTickLabel',tickLabel,'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
            colormap(hInhWidth,blueOnly);
        end

    filename = [dir,'VeffHalfWidth'];
    if strcmp(format,'fig')
        saveas(gcf,[filename,'.fig']);
    else
        print(gcf,[filename,'.',format],printDriver,dpi);
    end
end
