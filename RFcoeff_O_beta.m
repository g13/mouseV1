function coMat = RFcoeff_O_beta(lgnfile,nx,ny,nsig,save2file,tw,draw,format,ld,temper,reverse,OnOff,ie,antiphase,threads)
    if nargin < 15
        threads = 1;
        if nargin < 14
            antiphase = false;
            if nargin < 13
                ie = false;
                if nargin < 12
                    OnOff = false;
                    if nargin < 11
                        reverse = 1;
                        if nargin < 10
                            temper = 0.0;
                            if nargin < 9
                                ld = false;
                                if nargin < 8
                                    format = '';
                                    if nargin < 7
                                        draw = true;
                                        if nargin < 6
                                            tw = 1.0;
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
    addpath(genpath('./matlab_Utilities'));
    pPosition = [0, 0, 1280, 720];
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r100';
    end
    load(lgnfile);
    n = p.nv1e;
    m = p.nv1;
    maxd = max(p.enormDistance);
    mind = min(p.enormDistance);
    SepRFdist = maxd;
    %SepRFdist = max([0.64,max(p.enormDistance)]);
    set(0,'DefaultAxesFontSize',14);
    wE = p.enormDistance;
    wE(wE>=SepRFdist) = 1;
    wORF =wE(wE<SepRFdist);
    wE(wE<SepRFdist) = (wORF-mind)./SepRFdist;

    %SepRFdist = max([0.64,max(p.inormDistance)]);
    %wI = p.inormDistance;
    %wI(wI>=SepRFdist) = 1;
    %wORF =wI(wI<SepRFdist);
    %wI(wI<SepRFdist) = wORF./SepRFdist;

    switch reverse 
        case 1
            w = wE*tw;
        case 2
            w = ones(n,1)*tw;
        case -1 
            w = (1-wE)*tw;
        case -2
            w = ones(n,1)-tw;
    end
    oneMinusW = (1-w);  % oneMinusW * RFcorr + w * coeff
    wMat = ones(p.nv1,1)*w';
    oneMinusWmat = ones(p.nv1,1)*oneMinusW';

    if temper > 0 % for ORF, more overlap -> more blurness in dTheta gauge
        temperMat = ones(m,1) * (1-wE)' * temper;
        s = rng;
        temperMat = randn(m,n) .* temperMat * pi/180;
    else
        temperMat = zeros(m,n);
        s = -1;
    end

    clear wORF
    filename=['coMa-',num2str(nx),'x',num2str(ny),'-',lgnfile];
    global Aa Ab siga2 sigb2 X Y
    if ~ld
        AORF = false(p.nv1e,1);
        FWHM = 2*sqrt(2*log(2));
        Aa = 14.88 * (180/pi)^2;
        Ab = Aa*0.97;
        rc = 5.61; %deg
        rs = 16.98*rc/5.61; %deg
        siga = rc/sqrt(2);
        sigb = rs/sqrt(2);
        siga2 = siga^2;
        sigb2 = sigb^2;
        x = linspace(LGNpos(1,1)*180/pi - nsig*sigb,LGNpos(p.lgny*(p.lgnx-1)+1,1)*180/pi+nsig*sigb,nx);
        y = linspace(LGNpos(1,2)*180/pi - nsig*sigb,LGNpos(p.lgny,2)*180/pi+nsig*sigb,ny);
        ngrid = nx*ny;
        [X, Y] = meshgrid(x,y);
        rho = @(theta,a,b) (a.*b).^2./((b.*cos(theta)).^2+(a*sin(theta)).^2);
        pox = reshape(X,[ny*nx,1]);
        poy = reshape(Y,[ny*nx,1]);
        dpox = x(2)-x(1);
        dpoy = y(2)-y(1);
        ramp = sqrt(dpox^2 + dpoy^2);
        %ramp = 0;
        disp(['ramp size for grid involved: ',num2str(sqrt(ramp)),' degree']);
        figure;
        plot(LGNpos(:,1)*180/pi,LGNpos(:,2)*180/pi,'*k');
        xlim([min(x),max(x)]);
        ylim([min(y),max(y)]);
        Z = zeros(numel(X),m);
        tic;
        for i = 1:m
            for j = 1:length(nSubLGN{i})
                if nSubLGN{i}(j) > 0
                    Z(:,i) = Z(:,i) + multiLGNspatialKernel(LGNpos(v1Map{i,j},:)*180/pi,nSubLGN{i}(j),lgnStrength{i,j},1);
                else
                    Z(:,i) = Z(:,i) - multiLGNspatialKernel(LGNpos(-v1Map{i,j},:)*180/pi,-nSubLGN{i}(j),lgnStrength{i,j},1);
                end
            end
        end
        disp('RF grids ready');
        toc;
        coMat = ones(m,n);
        dThetaMat = ones(m,n);
        theta = [etheta;itheta];
        gtheta = zeros(m,1);
        pickTheta = theta>pi/2;
        gtheta(pickTheta) = theta(pickTheta) - pi/2;
        gtheta(~pickTheta) = theta(~pickTheta) + pi/2;
        RFcorrMat = ones(m,n);
        tic;
        nmZ = (Z-ones(ngrid,1)*mean(Z,1))./(ones(ngrid,1)*std(Z,1,1));
        disp('data normalized');
        toc;
        tic;
        parfor i=1:n
            if sum(nmZ(:,i)) == 0.0
                AORF(i) = true;
                RFcorrMat(:,i) = zeros(m,1);
            else
                jpick = [(i+1):m];
                        % [1, ngrid] x [ngird, m-i]
                RFcorr = (nmZ(:,i)'*nmZ(:,jpick))'/ngrid;
                RFcorr(isnan(RFcorr)) = 0;
                RFcorrMat(:,i) = [ones(i,1);RFcorr];
            end
            %if mod(i,round(n/10)) == 0
            %    disp(num2str(round(i/(n/10))));
            %end
        end
        disp('RFcoef ready');
        toc;
        if antiphase
            RFcorrMat((n+1):m) = -RFcorrMat((n+1):m);
        end
        for i=1:n
            RFcorrMat(i,(i+1):n) = RFcorrMat((i+1):n,i);
            RFcorrMat(i,i) = 1;
        end
        tic;
        thetaC = zeros(m,n);
        parfor i=1:n
            dTheta = abs(gtheta-gtheta(i));
            pickTheta = dTheta>pi/2;
            dTheta(pickTheta) = pi - dTheta(pickTheta);

            if temper > 0
                dTheta_tempered = dTheta + temperMat(:,i);
            else
                dTheta_tempered = dTheta;
            end
            
            dTheta_tempered = abs(dTheta_tempered);
            pickTheta = dTheta_tempered > pi/2;
            dTheta_tempered(pickTheta) = pi-dTheta_tempered(pickTheta);

            thetaC(:,i) = 1 - dTheta_tempered./(pi/4);
            dThetaMat(:,i) = dTheta;
        end
        disp('thetaC ready');
        toc;
        %coMat = RFcorrMat.*(ones(m,1)*w') + dThetaMat.*(ones(m,1)*oneMinusW');
        tic;
        coMat = RFcorrMat.*oneMinusWmat + thetaC.*wMat;
        disp('coMat finished');
        toc;
        nAORF = sum(AORF);
    else
        load([filename,'_more']);
        if temper > 0
            dThetaMat_tempered = dThetaMat + temperMat;
            dThetaMat_tempered = abs(dThetaMat_tempered);
            pickTheta = dThetaMat_tempered > pi/2;
            dThetaMat_tempered(pickTheta) = pi - dThetaMat_tempered(pickTheta);
        else
            dThetaMat_tempered = dThetaMat;
        end
        thetaC = 1-dThetaMat_tempered./(pi/4);
        coMat = RFcorrMat.*oneMinusWmat + thetaC.*wMat;
    end
    weights = figure;
    subplot(1,2,1)
    plot(p.enormDistance,oneMinusW,'o');
    title('1-w');
    subplot(1,2,2)
    plot(p.enormDistance,w,'o');
    title('w');
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(weights,['Weights-thetaC',num2str(tw*100),'%-',lgnfile,'.fig']);
        else
            print(weights,['Weights-thetaC',num2str(tw*100),'%-',lgnfile,'.',format],printDriver,dpi);
        end
    end

    disp(['# of Absolute ginger ale: ', num2str(nAORF)]);
    if save2file
        if ~ld
            save([filename,'_more'],'RFcorrMat','dThetaMat','AORF','nAORF','s');
        end
        save(filename,'coMat');
    end
    if draw
        ORF = p.typeE==1;
        SRF = p.typeE==2;
        nORF = sum(ORF); 
        nSRF = sum(SRF); 
        h = figure;
        dbin = 0.05;
        binranges = -1.0:dbin:1;
        nbins = length(binranges)-1;
        binranges(nbins+1) = 1.1;
        [binCountsE,indE] = histc(coMat(1:p.nv1e,:),binranges);       
        [binCountsI,indI] = histc(coMat(p.nv1e+(1:p.nv1i),:),binranges);       
        binCountsE = binCountsE(1:nbins,:);
        binCountsI = binCountsI(1:nbins,:);
        binranges = binranges(1:nbins)+dbin/2;

        subplot(2,4,1)
        hold on
        iE = indE(:,ORF);
        iI = indI(:,ORF);
        meanThetaE = zeros(nbins,1);
        stdThetaE = zeros(nbins,1);
        meanThetaI = zeros(nbins,1);
        stdThetaI = zeros(nbins,1);
        target = dThetaMat(1:p.nv1e,ORF);
        for i = 1:nbins
            meanThetaE(i) = mean(target(iE==i));
            stdThetaE(i) = std(target(iE==i));
        end
        target = dThetaMat(p.nv1e+(1:p.nv1i),ORF);
        for i = 1:nbins
            meanThetaI(i) = mean(target(iI==i));
            stdThetaI(i) = std(target(iI==i));
        end
        meanThetaE = meanThetaE.*180/pi;
        stdThetaE = stdThetaE .*180/pi;
        meanThetaI = meanThetaI.*180/pi;
        stdThetaI = stdThetaI .*180/pi;
        errorbar(binranges,meanThetaE,stdThetaE,'r');
        errorbar(binranges+dbin/5,meanThetaI,stdThetaI,'b');
        title('ORF')
        xlim([-1,1]);
        xlabel('Coefficient');
        ylabel('\Delta\theta')
        box on

        subplot(2,4,2)
        hold on
        iE = indE(:,SRF);
        iI = indI(:,SRF);
        meanThetaE = zeros(nbins,1);
        stdThetaE = zeros(nbins,1);
        meanThetaI = zeros(nbins,1);
        stdThetaI = zeros(nbins,1);
        target = dThetaMat(1:p.nv1e,SRF);
        for i = 1:nbins
            meanThetaE(i) = mean(target(iE==i));
            stdThetaE(i) = std(target(iE==i));
        end
        target = dThetaMat(p.nv1e+(1:p.nv1i),SRF);
        for i = 1:nbins
            meanThetaI(i) = mean(target(iI==i));
            stdThetaI(i) = std(target(iI==i));
        end
        meanThetaE = meanThetaE.*180/pi;
        stdThetaE = stdThetaE .*180/pi;
        meanThetaI = meanThetaI.*180/pi;
        stdThetaI = stdThetaI .*180/pi;
        errorbar(binranges,meanThetaE,stdThetaE,'r');
        errorbar(binranges+dbin/5,meanThetaI,stdThetaI,'b');
        xlim([-1,1]);
        xlabel('Coefficient');
        ylabel('\Delta\theta')
        title('SRF');
        box on

        subplot(2,2,3)
        title('Dist. of Coeffs');
        bar(binranges,[sum(binCountsE(:,ORF),2),sum(binCountsE(:,SRF),2),sum(binCountsI(:,ORF),2),sum(binCountsI(:,SRF),2)]);
        legend({'ORF_E','SRF_E','ORF_I','SRF_I'});
        xlabel('Coefficient');
        ylabel('# Cells');

        dbin = 0.05;
        binranges = -1.0:dbin:1;
        nbins = length(binranges)-1;
        binranges(nbins+1) = 1.1;
        [binCountsE,indE] = histc(RFcorrMat(1:p.nv1e,:),binranges);       
        [binCountsI,indI] = histc(RFcorrMat(p.nv1e+(1:p.nv1i),:),binranges);       
        binCountsE = binCountsE(1:nbins,:);
        binCountsI = binCountsI(1:nbins,:);
        binranges = binranges(1:nbins)+dbin/2;

        subplot(2,4,3)
        hold on
        iE = indE(:,ORF);
        iI = indI(:,ORF);
        meanThetaE = zeros(nbins,1);
        stdThetaE = zeros(nbins,1);
        meanThetaI = zeros(nbins,1);
        stdThetaI = zeros(nbins,1);
        target = dThetaMat(1:p.nv1e,ORF);
        for i = 1:nbins
            meanThetaE(i) = mean(target(iE==i));
            stdThetaE(i) = std(target(iE==i));
        end
        target = dThetaMat(p.nv1e+(1:p.nv1i),ORF);
        for i = 1:nbins
            meanThetaI(i) = mean(target(iI==i));
            stdThetaI(i) = std(target(iI==i));
        end
        meanThetaE = meanThetaE.*180/pi;
        stdThetaE = stdThetaE .*180/pi;
        meanThetaI = meanThetaI.*180/pi;
        stdThetaI = stdThetaI .*180/pi;
        errorbar(binranges,meanThetaE,stdThetaE,'r');
        errorbar(binranges+dbin/5,meanThetaI,stdThetaI,'b');
        title('ORF');
        xlim([-1,1]);
        xlabel('RFcorr');
        ylabel('\Delta\theta')
        box on

        subplot(2,4,4)
        hold on
        iE = indE(:,SRF);
        iI = indI(:,SRF);
        meanThetaE = zeros(nbins,1);
        stdThetaE = zeros(nbins,1);
        meanThetaI = zeros(nbins,1);
        stdThetaI = zeros(nbins,1);
        target = dThetaMat(1:p.nv1e,SRF);
        for i = 1:nbins
            meanThetaE(i) = mean(target(iE==i));
            stdThetaE(i) = std(target(iE==i));
        end
        target = dThetaMat(p.nv1e+(1:p.nv1i),SRF);
        for i = 1:nbins
            meanThetaI(i) = mean(target(iI==i));
            stdThetaI(i) = std(target(iI==i));
        end
        meanThetaE = meanThetaE.*180/pi;
        stdThetaE = stdThetaE .*180/pi;
        meanThetaI = meanThetaI.*180/pi;
        stdThetaI = stdThetaI .*180/pi;
        errorbar(binranges,meanThetaE,stdThetaE,'r');
        errorbar(binranges+dbin/5,meanThetaI,stdThetaI,'b');
        title('SRF');
        xlim([-1,1]);
        xlabel('RFcorr');
        ylabel('\Delta\theta')
        box on

        subplot(2,2,4)
        title('Dist. of RFcorr');
        bar(binranges,[sum(binCountsE(:,ORF),2),sum(binCountsE(:,SRF),2),sum(binCountsI(:,ORF),2),sum(binCountsI(:,SRF),2)]);
        legend({'ORF_E','SRF_E','ORF_I','SRF_I'});
        xlabel('RFcorr');
        ylabel('# Cells');

        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(h,['Pop_O_beta_RF-dTheta',num2str(tw*100),'%-',lgnfile,'.fig']);
            else
                print(h,['Pop_O_beta_RF-dTheta',num2str(tw*100),'%-',lgnfile,'.',format],printDriver,dpi);
            end
        end
%% on dTheta
        thetaC = 1-dThetaMat./(pi/4);
        hHeat = figure;
        dCV = 0.05;
        ctrs = cell(2,1);
        ctrs{1} = linspace(-1,1,7);
        ctrs{2} = -1:dCV:1;
        lctrsx = length(ctrs{1});
        lctrsy = length(ctrs{2});
        dTickX = 1/lctrsx;
        dTickY = 0.1;
        tickPosX = linspace(0.5,lctrsx-1+0.5,lctrsx);
        tickPosY = 0.5:lctrsy*dTickY:lctrsy+0.5;
        tickLabelX = num2str((linspace(0,90,length(tickPosX)))');
        tickLabelY = flipud(num2str((linspace(-1,1,length(tickPosY)))'));

        hExc = subplot(2,2,1);
        RFpair = [reshape(thetaC(1:p.nv1e,ORF),[p.nv1e*nORF,1]), reshape(coMat(1:p.nv1e,ORF),[p.nv1e*nORF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('Coefficient');
        title('Exc -> ORF');
        colormap(hExc,redOnly);

        hInh = subplot(2,2,2);
        RFpair = [reshape(thetaC(p.nv1e+(1:p.nv1i),ORF),[p.nv1i*nORF,1]), reshape(coMat(p.nv1e+(1:p.nv1i),ORF),[p.nv1i*nORF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('Coefficient');
        title('Inh -> ORF');
        colormap(hInh,blueOnly);

        hExc = subplot(2,2,3);
        RFpair = [reshape(thetaC(1:p.nv1e,SRF),[p.nv1e*nSRF,1]), reshape(coMat(1:p.nv1e,SRF),[p.nv1e*nSRF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx,1:lctrsy-1);
        %denRFpair = denRFpair(1:(lctrsx-1,
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('Coefficient');
        title('Inh -> SRF');
        colormap(hExc,redOnly);

        hInh = subplot(2,2,4);
        RFpair = [reshape(thetaC(p.nv1e+(1:p.nv1i),SRF),[p.nv1i*nSRF,1]), reshape(coMat(p.nv1e+(1:p.nv1i),SRF),[p.nv1i*nSRF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('Coefficient');
        title('Inh -> SRF');
        colormap(hInh,blueOnly);
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(hHeat,['Pop_O_beta_RF-Coeff_vs_dTheta',num2str(tw*100),'%-',lgnfile,'.fig']);
            else
                print(hHeat,['Pop_O_beta_RF-Coeff_vs_dTheta',num2str(tw*100),'%-',lgnfile,'.',format],printDriver,dpi);
            end
        end
        hHeat = figure;
        dCV = 0.05;
        ctrs = cell(2,1);
        ctrs{1} = linspace(-1,1,7);
        ctrs{2} = -1:dCV:1;
        lctrsx = length(ctrs{1});
        lctrsy = length(ctrs{2});
        dTickX = 2/lctrsx;
        dTickY = 0.1;
        tickPosX = linspace(0.5,lctrsx-1+0.5,lctrsx);
        tickPosY = 0.5:lctrsy*dTickY:lctrsy+0.5;
        tickLabelX = num2str((linspace(0,90,length(tickPosX)))');
        tickLabelY = flipud(num2str((linspace(-1,1,length(tickPosY)))'));

        hExc = subplot(2,2,1);
        RFpair = [reshape(thetaC(1:p.nv1e,ORF),[p.nv1e*nORF,1]), reshape(RFcorrMat(1:p.nv1e,ORF),[p.nv1e*nORF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('RFcorr');
        title('Exc -> ORF');
        colormap(hExc,redOnly);

        hInh = subplot(2,2,2);
        RFpair = [reshape(thetaC(p.nv1e+(1:p.nv1i),ORF),[p.nv1i*nORF,1]), reshape(RFcorrMat(p.nv1e+(1:p.nv1i),ORF),[p.nv1i*nORF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('RFcorr');
        title('Inh -> ORF');
        colormap(hInh,blueOnly);

        hExc = subplot(2,2,3);
        RFpair = [reshape(thetaC(1:p.nv1e,SRF),[p.nv1e*nSRF,1]), reshape(RFcorrMat(1:p.nv1e,SRF),[p.nv1e*nSRF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('RFcorr');
        title('Inh -> SRF');
        colormap(hExc,redOnly);

        hInh = subplot(2,2,4);
        RFpair = [reshape(thetaC(p.nv1e+(1:p.nv1i),SRF),[p.nv1i*nSRF,1]), reshape(RFcorrMat(p.nv1e+(1:p.nv1i),SRF),[p.nv1i*nSRF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx],[1,lctrsy-1],denRFpair');
        hold on
        plot(linspace(lctrsx+0.5,0.5,5),linspace(0.5,lctrsy+0.5,5),'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('RFcorr');
        title('Inh -> SRF');
        colormap(hInh,blueOnly);
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(hHeat,['Pop_O_beta_RFcorr_vs_dTheta',num2str(tw*100),'%-',lgnfile,'.fig']);
            else
                print(hHeat,['Pop_O_beta_RFcorr_vs_dTheta',num2str(tw*100),'%-',lgnfile,'.',format],printDriver,dpi);
            end
        end
%% heat Map
        hHeat = figure;
        dCV = 0.02;
        ctrs = cell(2,1);
        ctrs{1} = -1:dCV:1-dCV;
        ctrs{2} = ctrs{1};
        lctrsx = length(ctrs{1});
        lctrsy = length(ctrs{2});
        dTick = 0.1;
        tickPosY = 0.5:lctrsx*dTick:lctrsx+0.5;
        tickPosX = 0.5:lctrsy*dTick:lctrsy+0.5;
        tickLabel = num2str((linspace(-1,1,length(tickPosX)))');

        hExc = subplot(2,2,1);
        RFpair = [reshape(RFcorrMat(1:p.nv1e,ORF),[p.nv1e*nORF,1]), reshape(coMat(1:p.nv1e,ORF),[p.nv1e*nORF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx-1],[lctrsy-1,1],denRFpair');
        hold on
        plot(lctrsx+0.5:-1:0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
        xlabel('RFcorr');
        ylabel('Coefficient');
        title('Exc -> ORF');
        colormap(hExc,redOnly);

        hInh = subplot(2,2,2);
        RFpair = [reshape(RFcorrMat(p.nv1e+(1:p.nv1i),ORF),[p.nv1i*nORF,1]), reshape(coMat(p.nv1e+(1:p.nv1i),ORF),[p.nv1i*nORF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx-1],[lctrsy-1,1],denRFpair');
        hold on
        plot(lctrsx+0.5:-1:0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
        xlabel('RFcorr');
        ylabel('Coefficient');
        title('Inh -> ORF');
        colormap(hInh,blueOnly);

        hExc = subplot(2,2,3);
        RFpair = [reshape(RFcorrMat(1:p.nv1e,SRF),[p.nv1e*nSRF,1]), reshape(coMat(1:p.nv1e,SRF),[p.nv1e*nSRF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx-1],[lctrsy-1,1],denRFpair');
        hold on
        plot(lctrsx+0.5:-1:0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
        xlabel('RFcorr');
        ylabel('Coefficient');
        title('Inh -> SRF');
        colormap(hExc,redOnly);

        hInh = subplot(2,2,4);
        RFpair = [reshape(RFcorrMat(p.nv1e+(1:p.nv1i),SRF),[p.nv1i*nSRF,1]), reshape(coMat(p.nv1e+(1:p.nv1i),SRF),[p.nv1i*nSRF,1])];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx-1],[lctrsy-1,1],denRFpair');
        hold on
        plot(lctrsx+0.5:-1:0.5, 0.5:lctrsy+0.5,'-.k','LineWidth',2);
        
        set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
        xlabel('RFcorr');
        ylabel('Coefficient');
        title('Inh -> SRF');
        colormap(hInh,blueOnly);
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(hHeat,['Pop_O_beta_RFvsCoeff-dTheta',num2str(tw*100),'%-',lgnfile,'.fig']);
            else
                print(hHeat,['Pop_O_beta_RFvsCoeff-dTheta',num2str(tw*100),'%-',lgnfile,'.',format],printDriver,dpi);
            end
        end
        if nAORF > 0
        hORF = figure;
        subplot(2,1,1)
        hold on
        sel = AORF;
        nsel = sum(sel);
        binranges = -1:0.1:1;
        nbins = length(binranges);
        binranges(nbins) = 1.1;
        sbC = zeros(nbins-1,4);

        tmp = coMat(1:p.nv1e,sel);
        bC = histc(tmp,binranges);
        sbC(:,1) = sum(bC(1:nbins-1,:),2);

        tmp = RFcorrMat(1:p.nv1e,sel);
        bC = histc(tmp,binranges);
        sbC(:,2) = sum(bC(1:nbins-1,:),2);

        tmp = coMat(p.nv1e+(1:p.nv1i),sel);
        bC = histc(tmp,binranges);
        sbC(:,3) = sum(bC(1:nbins-1,:),2);

        tmp = RFcorrMat(p.nv1e+(1:p.nv1i),sel);
        bC = histc(tmp,binranges);
        sbC(:,4) = sum(bC(1:nbins-1,:),2);
        plot(binranges(1:nbins-1),sbC,'*','LineStyle','-');

        legend({'coMatE','RFcorrE','coMatI','RFcorrI'});

        hExc = subplot(2,2,3);
        dCoeff = 0.05;
        ctrs = cell(2,1);
        ctrs{1} = linspace(0,pi/2,7);
        ctrs{2} = -1:dCoeff:1;
        lctrsx = length(ctrs{1});
        lctrsy = length(ctrs{2});
        ctrs{2}(end) = 1+dCoeff;
        dTickX = 1/lctrsx;
        TickY = 0.1;
        tickPosX = linspace(0.5,lctrsx-1+0.5,lctrsx);
        tickPosY = 0.5:lctrsy*dTickY:lctrsy+0.5;
        tickLabelX = num2str((linspace(0,90,length(tickPosX)))');
        tickLabelY = flipud(num2str((linspace(-1,1,length(tickPosY)))'));

        tmp = coMat(1:p.nv1e,sel);
        tmp = reshape(tmp,[p.nv1e*nsel,1]);
        tmpTheta = dThetaMat(1:p.nv1e,sel);
        tmpTheta = reshape(tmpTheta,[p.nv1e*nsel,1]);
        RFpair = [tmpTheta,tmp];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx-1],[lctrsy-1,1],denRFpair');
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('Coeff');
        title('exc');
        colormap(hExc,redOnly);
        hInh = subplot(2,2,4);

        tmp = coMat(p.nv1e+(1:p.nv1i),sel);
        tmp = reshape(tmp,[p.nv1i*nsel,1]);
        tmpTheta = dThetaMat(p.nv1e+(1:p.nv1i),sel);
        tmpTheta = reshape(tmpTheta,[p.nv1i*nsel,1]);
        RFpair = [tmpTheta,tmp];
        denRFpair = hist3(RFpair,ctrs);
        denRFpair = denRFpair(1:lctrsx-1,1:lctrsy-1);
        maxDen = max(max(denRFpair));
        denRFpair = denRFpair/maxDen;
        imagesc([1,lctrsx-1],[lctrsy-1,1],denRFpair');
        set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
        xlabel('\Delta\theta');
        ylabel('Coeff');
        title('inh');
        colormap(hInh,blueOnly);
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                saveas(hORF,['Pop_O_beta_ORF-dTheta',num2str(tw*100),'%-',lgnfile,'.fig']);
            else
                print(hORF,['Pop_O_beta_ORF-dTheta',num2str(tw*100),'%-',lgnfile,'.',format],printDriver,dpi);
            end
        end
        end
    end
end

function Z = multiLGNspatialKernel(pos,n,s,scale)
    global Aa Ab siga2 sigb2 X Y
    Z = zeros(size(X));
    for i = 1:n
        Z = Z + s(i)*(Aa/(siga2*2*pi) * exp(-((X-pos(i,1)).^2+(Y-pos(i,2)).^2)/siga2/2)...
            - Ab/(sigb2*2*pi) * exp(-((X-pos(i,1)).^2+(Y-pos(i,2)).^2)/sigb2/2));
    end
    Z = Z/scale;
    Z = reshape(Z,[numel(Z),1]);
%     assert(sum(Z)>0);
end
