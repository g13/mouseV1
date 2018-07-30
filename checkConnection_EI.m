function checkConnection_EI(meefile,lgnfile,coMatfile,checklist,nn,format,draw,ntheta,Inh,spread,theme)
    global Aa Ab siga2 sigb2 x y X Y LGNpos v1Map lgnStrength nSubLGN p etheta itheta rfx rfy ttt FWHM plotPr
    pPosition = [0, 0, 1280, 720];
    plotPr = 5;
    if nargin < 10
        spread = false;
        if nargin < 9
            Inh = false;
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
    set(0,'DefaultAxesFontSize',14);
    load(lgnfile);
    load([coMatfile,'_more'],'RFcorrMat','dThetaMat');
    load(coMatfile,'coMat');
    coMatE = coMat(1:p.nv1e,1:p.nv1e);
    RFcorrMatE = RFcorrMat(1:p.nv1e,1:p.nv1e);
    dThetaE = dThetaMat(1:p.nv1e,1:p.nv1e);
    if Inh
        coMatI = coMat(p.nv1e+(1:p.nv1i),1:p.nv1e);
        dThetaI = dThetaMat(p.nv1e+(1:p.nv1i),1:p.nv1e);
        RFcorrMatI = RFcorrMat(p.nv1e +(1:p.nv1i),1:p.nv1e);
    end
    clear RFcorrMatI; 
    clear dThetaMat;
    load(meefile);
    meifile = meefile;
    meifile(3) = 'i';
    load(meifile);
    load('logNormalProfile.mat');
    FontSize = 11;
    set(0,'DefaultAxesFontSize',FontSize);
    checklist = reshape(checklist,[numel(checklist),1]);
 
    rfx = @(x0,a,b,theta,t,ra,rb) x0 + (ra.*a.*cos(t).*cos(theta)-rb.*b.*sin(t).*sin(theta));
    rfy = @(y0,a,b,theta,t,ra,rb) y0 + (ra.*a.*cos(t).*sin(theta)+rb.*b.*sin(t).*cos(theta));
    ttt = 0:0.1:2*pi;
    FWHM = 2*sqrt(2*log(2));

    Aa = 14.88 * (180/pi)^2;
    Ab = Aa*0.97;
    rc = 5.61; %deg
    rs = 16.98*rc/5.61; %deg
    siga = rc/sqrt(2);
    sigb = rs/sqrt(2);
    siga2 = siga^2;
    sigb2 = sigb^2;
    nsig = 3;
    nx = 30;
    ny = 30;
    x = linspace(LGNpos(1,1)*180/pi - nsig*sigb,LGNpos(p.lgny*(p.lgnx-1)+1,1)*180/pi+nsig*sigb,nx);
    y = linspace(LGNpos(1,2)*180/pi - nsig*sigb,LGNpos(p.lgny,2)*180/pi+nsig*sigb,ny);
    [X, Y] = meshgrid(x,y);
    nlist = length(checklist);
    mcoeff = zeros(nlist,1);
    scoeff = mcoeff;
    dbin = 0.05;
    dtheta = pi/ntheta;
    bound = sqrt(2);

    hCossell = figure;
    pickedNeighbor = mee > 0;
    RFflat = RFcorrMatE;
    RFflat_connected = RFcorrMatE(pickedNeighbor);
    excStr = zeros(p.nv1e);
    if spread
        for i = 1:p.nv1e
            ineighbor = pickedNeighbor(:,i);
            excStr(ineighbor,i) = profiles(mee(ineighbor,i),i);
        end
    else
        excStr(pickedNeighbor) = profiles(mee(pickedNeighbor));
    end
    excStr_connected = excStr(pickedNeighbor);
    
    [RFflatC_sorted, ind] = sort(RFflat_connected);
    excStrC_sorted = excStr_connected(ind);
    excWeight_count = cumsum(excStr_sorted);
    connected_count = cumsum(ones(length(excStr_connected),1));
    dbin = 0.01;
    binranges = -1.0:dbin:1;
    nbins = length(binranges);
    binranges(nbins+1) = 1+dbin;
    [RFflat_count, ~] = histc(RFflat,binranges)
    binranges = binranges(1:nbins)+dbin/2;
    subplot(1,2,1)
    bar(binranges,RFflat_count);
    subplot(1,2,2)
    hold on
    plot(RFlatC_sorted,excWeight_count);
    plot(RFlatC_sorted,connected_count);
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(h,['Cossell-',coMatfile,'-',theme,'.',format],printDriver,dpi);
        saveas(h,['Cossell-',coMatfile,'-',theme,'.fig']);
    end

    if nlist==p.nv1e
        h = figure;
        subplot(1,2,1)
        pickedNeighbor = mee>0;
        binranges = -1.0:dbin:1;
        nbins = length(binranges)-1;
        binranges(nbins+1) = 1.1;
        [binCounts0,ind] = histc(coMatE,binranges);       
        binCounts0 = binCounts0(1:nbins,1:p.nv1e);
        binCounts0 = binCounts0./(ones(nbins,1)*sum(binCounts0,1))*100;
        binCounts0 = mean(binCounts0,2);
        
        coMatE(~pickedNeighbor) = binranges(1)-1;
        [binCounts,~] = histc(coMatE,binranges);
        binCounts = binCounts(1:nbins,1:p.nv1e);
        binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
        binCounts = mean(binCounts,2);
        
        EPSPslice = zeros(nbins,1);
        for j = 1:nbins
            if spread
                for i = 1:p.nv1e
                    ineighbor = mee(:,i)>0;
                    EPSPslice(j) = EPSPslice(j) + sum(profiles(mee((ind(:,i)==j)&ineighbor,i),i));
                end
                EPSPslice(j) = EPSPslice(j)/p.nv1e;
            else
                EPSPslice(j) = sum(profiles(mee((ind==j)&pickedNeighbor)))/p.nv1e;
            end
        end
        EPSPslice = EPSPslice/sum(EPSPslice)*100;
        binranges = binranges(1:nbins)+dbin/2;
        plot(binranges,[cumsum(EPSPslice),cumsum(binCounts),binCounts0./sum(binCounts0)*100]);
        ylabel('%');
        xlabel('coMat');
        legend({'Cort. Exc. CDF.','Connected CDF.','PDF'},'Location','NorthWest');
        ylim([0,100]);

        subplot(1,2,2)
        binranges = -1.0:dbin:1;
        nbins = length(binranges)-1;
        binranges(nbins+1) = 1.1;
        [binCounts0,ind] = histc(RFcorrMatE,binranges);       
        binCounts0 = binCounts0(1:nbins,1:p.nv1e);
        binCounts0 = binCounts0./(ones(nbins,1)*sum(binCounts0,1))*100;
        binCounts0 = mean(binCounts0,2);
        
        RFcorrMatE(~pickedNeighbor) = binranges(1)-1;
        [binCounts,~] = histc(RFcorrMatE,binranges);
        binCounts = binCounts(1:nbins,1:p.nv1e);
        binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
        binCounts = mean(binCounts,2);
        
        EPSPslice = zeros(nbins,1);
        for j = 1:nbins
            if spread
                for i = 1:p.nv1e
                    ineighbor = pickedNeighbor(:,i);
                    EPSPslice(j) = EPSPslice(j) + sum(profiles(mee((ind(:,i)==j)&ineighbor,i),i));
                end
                EPSPslice(j) = EPSPslice(j)/p.nv1e;
            else
                EPSPslice(j) = sum(profiles(mee((ind==j)&pickedNeighbor)))/p.nv1e;
            end
        end
        EPSPslice = EPSPslice/sum(EPSPslice)*100;
        binranges = binranges(1:nbins)+dbin/2;
        plot(binranges,[cumsum(EPSPslice),cumsum(binCounts),binCounts0./sum(binCounts0)*100]);
        ylabel('%');
        xlabel('RF correlation');
        legend({'Cort. Exc. CDF.','Connected CDF.','PDF'},'Location','NorthWest');
        ylim([0,100]);
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(h,['Pop_RFsimilarity-',coMatfile,'-',theme,'.',format],printDriver,dpi);
            saveas(h,['Pop_RFsimilarity-',coMatfile,'-',theme,'.fig']);
        end
        h = figure;
    %% theta
            
        binranges = 0.0:dtheta:(pi/2);
        nbins = length(binranges)-1;
        binranges(nbins+1) = pi/2+0.1;
        tt = {'ORF','SRF'};
        for type = 1:2
            subplot(2,2,type);
            sel = p.typeE==type;
            meeT = mee(:,sel);
            pickedNeighborT = pickedNeighbor(:,sel);
            nSel = sum(sel);
            dThetaMatC = dThetaE(:,sel);
            if spread
                profilesT = profiles(:,sel);
            end
            dThetaMatC(~pickedNeighborT) = binranges(1)-1;
            [binCounts,ind] = histc(dThetaMatC,binranges);
            binCounts = binCounts(1:nbins,:);
            binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
            mbinCounts = mean(binCounts,2);
            mbinCounts = mbinCounts./sum(mbinCounts)*100;
            stdbinCounts = std(binCounts,1,2);

            EPSPslice = zeros(nbins,1);
            mEPSPslice = zeros(nbins,1);
            for j = 1:nbins
                if spread
                for i = 1:nSel
                    ineighbor = pickedNeighborT(:,i);
                    tmp = profilesT(meeT((ind(:,i)==j)&ineighbor,i),i);
                    EPSPslice(j) = EPSPslice(j) + sum(tmp);
                    mEPSPslice(j) = mEPSPslice(j) + mean(tmp);
                end
                mEPSPslice(j) = mEPSPslice(j)/nSel;
                else
                    tmp = profiles(meeT((ind==j)&pickedNeighborT));
                    EPSPslice(j) = sum(tmp);
                    mEPSPslice(j) = mean(tmp);
                end
            end
            EPSPslice = EPSPslice./sum(EPSPslice)*100;
            xxxx = (binranges(1:nbins)+dtheta/2)*180/pi;
            [hAx,h1,h2] = plotyy(xxxx,mEPSPslice,xxxx,[mbinCounts,EPSPslice]);
            h1.Marker = '.';
            h2(1).LineStyle = ':';
            h2(1).Marker = 'o';
            h2(2).Marker = 's';
            ylabel(hAx(1),'EPSP (mV)');
            ylabel(hAx(2),'%');
            ylim(hAx(1),[0,inf]);
            ylim(hAx(2),[0,inf]);
            xlabel('\Delta\theta');
            legend({'avg EPSP','% Connection','% Total EPSP'});
            title(tt(type));
        end
        subplot(2,1,2);
        pickedNeighbor = mei>0;
        dThetaMatC = dThetaI;
        dThetaMatC(~pickedNeighbor) = binranges(1)-1;
        binCounts = histc(dThetaMatC,binranges);
        binCounts = binCounts(1:nbins,:);
        binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
        mbinCounts = mean(binCounts,2);
        mbinCounts = mbinCounts./sum(mbinCounts)*100;
        stdbinCounts = std(binCounts,1,2);

        xxxx = (binranges(1:nbins)+dtheta/2)*180/pi;
        errorbar(xxxx,mbinCounts,stdbinCounts);
        ylim([0,inf]);
        ylabel('%');
        xlabel('\Delta\theta');
        legend('% Connection');
        title('inh');

        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(h,['Pop_Theta-',coMatfile,'-',theme,'.',format],printDriver,dpi);
            saveas(h,['Pop_Theta-',coMatfile,'-',theme,'.fig']);
        end

       return
    end
    if bound> 0
        e.raxn = 100;
        e.rden = 75;
        iii.raxn = 100;
        iii.rden = 75;
        elx = (p.ev1x) * p.dE;
        ely = (p.ev1y) * p.dE;
        pex = ones(p.ev1y,1)*((1:p.ev1x)-0.5).*p.dE;
        pex = reshape(pex,[p.nv1e,1]); % id goes along y first (column major)
        pey = ((1:p.ev1y)-0.5)'*ones(1,p.ev1x).*p.dE;
        pey = reshape(pey,[p.nv1e,1]);

        src.x = pex;
        src.y = pey;
        src.lx = elx;
        src.ly = ely;
        src.raxn = e.raxn;
        src.n = p.nv1e;
        tar.x = pex;
        tar.y = pey;
        tar.n = p.nv1e;
        tar.rden = e.rden;
        tar.bound = bound;
        cR = checkRegion(checklist,src,tar);
    else
        cR = true(p.nv1e,nlist);
    end
    for k=1:nlist
        i = checklist(k);
        pickedNeighbor = mee(:,i) > 0;
        pickedNeighbor(i) = false;

        coeffs = coMatE(:,i);
        RFcoeffs = RFcorrMatE(:,i);
        dThetas = dThetaE(:,i);
        [~,bigShots] = sort(mee(:,i),'descend');
        bigShots = bigShots(1:nn);
        if draw
            h = figure;
            subplot(ceil((nn-1)/plotPr)+1,3,1);
            
            binranges = -1.0:dbin:1;
            nbins = length(binranges)-1;
            binranges(nbins+1) = 1.1;
            binCounts = histc(RFcoeffs(pickedNeighbor),binranges);
            binCounts = binCounts(1:nbins);

            binCounts0 = histc(RFcoeffs(cR(:,k)),binranges);
            binCounts0 = binCounts0(1:nbins);

            if spread
                EPSPs = profiles(mee(pickedNeighbor,i),i);
            else
                EPSPs = profiles(mee(pickedNeighbor,i));
            end
            [Ax,h1,h2] = plotyy(RFcoeffs(pickedNeighbor),EPSPs,binranges(1:nbins)+dbin/2,[100*binCounts./sum(binCounts),100*binCounts0./sum(binCounts0)]);
            h1.LineStyle = 'none';
            h2(1).LineStyle = '--';
            h2(2).LineStyle = ':';
            h1.Marker = '.';
            h2(1).Marker = 'o';
            h2(2).Marker = 'o';
            ylabel(Ax(1),'EPSP (mV)');
            ylabel(Ax(2),'%');
            xlabel('RFcorr');
            legend({'EPSP','Pre.Neurons','Neuron Pool'},'Location','NorthWest');
            n50 = sum(binCounts((binranges(1:nbins)+dbin/2)>=0.5));
            title({[num2str(n50),' RFcorr >= 0.5 ',num2str(n50/sum(binCounts)*100,'%2.1f'),'%'],['EPSP =', num2str(mean(EPSPs)),'\pm',num2str(std(EPSPs))]});
            
            subplot(ceil((nn-1)/plotPr)+1,3,2);

            binranges = -1.0:dbin:1;
            nbins = length(binranges)-1;
            binranges(nbins+1) = 1.1;
            binCounts = histc(coeffs(pickedNeighbor),binranges);
            binCounts = binCounts(1:nbins);

            binCounts0 = histc(coeffs(cR(:,k)),binranges);
            binCounts0 = binCounts0(1:nbins);
            [~,ind] = histc(coeffs,binranges);
            
            binranges = binranges(1:nbins)+dbin/2;
            norm_binCounts = binCounts/sum(binCounts)*100;
            norm_binCounts0 = binCounts0/sum(binCounts0)*100;
            EPSPslice = zeros(nbins,1);
            for j = 1:nbins
                if spread
                    EPSPslice(j) = sum(profiles(mee((ind==j)&pickedNeighbor,i),i));
                else
                    EPSPslice(j) = sum(profiles(mee((ind==j)&pickedNeighbor,i)));
                end
            end
            EPSPslice = EPSPslice/sum(EPSPslice)*100;
            plotyy(binranges,[EPSPslice,norm_binCounts,norm_binCounts0],binranges,binCounts*100./binCounts0);
            ylabel('%');
            legend({'Cort. Exc.','Pre. Neuron','All Neurons','Prob'},'Location','NorthWest');
            ylim([0,inf]);
            xlabel('Coeff');
            check([i;bigShots],h,coeffs);

        %% theta
            subplot(ceil((nn-1)/plotPr)+1,3,9);
            binranges = 0.0:dtheta:(pi/2);
            nbins = length(binranges)-1;
            binranges(nbins+1) = pi/2+0.1;
            dThetas(~pickedNeighbor) = binranges(1)-1;
            [binCounts,ind] = histc(dThetas,binranges);
            binCounts = binCounts(1:nbins);
            mbinCounts = binCounts/sum(binCounts)*100;
            %EPSPs = profiles(mee(pickedNeighbor,i));
            EPSPslice = zeros(nbins,1);
            mEPSPslice = zeros(nbins,1);
            for j = 1:nbins
                %EPSPslice(j) = sum(profiles(mee((ind==j)&pickedNeighbor,i)));
                if spread
                    tmp = profiles(mee((ind==j)&pickedNeighbor,i),i);
                else
                    tmp = profiles(mee((ind==j)&pickedNeighbor,i));
                end
                mEPSPslice(j) = mean(tmp);
                EPSPslice(j) = sum(tmp);
            end
            EPSPslice = EPSPslice./sum(EPSPslice)*100;
            xxxx = (binranges(1:nbins)+dtheta/2)*180/pi;
            [hAx,h1,h2] = plotyy(xxxx,mEPSPslice,xxxx,[mbinCounts,EPSPslice]);
            h1.Marker = '.';
            h2(1).LineStyle = ':';
            h2(1).Marker = 'o';
            h2(2).Marker = 's';
            ylabel(hAx(1),'EPSP (mV)');
            ylabel(hAx(2),'%');
            legend({'avg EPSP','% Connection','% Total EPSP'});
            xlabel('\Delta\theta');
            if ~isempty(format)
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                print(h,['EE_RFsimilarity-',num2str(i),'-',coMatfile,'-',theme,'.',format],printDriver,dpi);
            end
        end
        mcoeff(i) = mean(coeffs(bigShots));
        scoeff(i) = std(coeffs(bigShots));
    end
    fp = 3;
    disp(['mean coeff ',num2str(mean(mcoeff),fp),'+-',num2str(std(mcoeff),fp)]);
    disp(['std of coeff ',num2str(mean(scoeff),fp),'+-',num2str(std(scoeff),fp)]);
    if Inh
        if bound > 0
            ilx = (p.iv1x) * p.dI;
            ily = (p.iv1y) * p.dI;
            pix = ones(p.iv1y,1) * (((1:p.iv1x)-0.5).*p.dI);
            pix = reshape(pix,[p.nv1i,1]);
            piy = (((1:p.iv1y)-0.5)'.*p.dI) * ones(1,p.iv1x);
            piy = reshape(piy,[p.nv1i,1]);

            src.x = pix;
            src.y = piy;
            src.n = p.nv1i;
            src.lx = ilx;
            src.ly = ily;
            src.raxn = iii.raxn;
            cR = checkRegion(checklist,src,tar);
        else
            cR = true(p.nv1e,nlist);
        end
        for k=1:nlist
            i = checklist(k);
            pickedNeighbor = mei(:,i) > 0;
            coeffs = coMatI(:,i);
            dThetas = dThetaI(:,i);
            if sum(pickedNeighbor) ~= sum(mei(:,i))
                [~,bigShots] = sort(mei(:,i),'descend');
                bigShots = bigShots(1:nn);
            else
                [~,bigShots] = sort(coeffs,'descend');
                sPN = pickedNeighbor(bigShots);
                bigShots = bigShots(sPN);
                bigShots = bigShots(1:nn);
            end
            if draw
                h = figure;
                subplot(ceil((nn-1)/plotPr)+1,3,1);
                plot(coeffs(pickedNeighbor),profiles(mei(pickedNeighbor,i)),'.');
                ylabel('IPSP (mV)');
                xlabel('RFcorr');
                title({['largest RFcorr:',num2str(max(coeffs))],...
                    ['largest Picked',num2str(max(coeffs(pickedNeighbor)))]});
                
                subplot(ceil((nn-1)/plotPr)+1,3,2);
                binranges = -1.0:dbin:1;
                nbins = length(binranges)-1;
                binranges(nbins+1) = 1.1;
                binCounts0 = histc(coeffs(cR(:,k)),binranges);
                binCounts0 = binCounts0(1:nbins);
                [~,ind] = histc(coeffs,binranges);
                binCounts = histc(coeffs(pickedNeighbor),binranges);
                binCounts = binCounts(1:nbins);
                binranges = binranges(1:nbins)+dbin/2;
                norm_binCounts = binCounts/sum(binCounts)*100;
                norm_binCounts0 = binCounts0/sum(binCounts0)*100;
                IPSPslice = zeros(nbins,1);
                for j = 1:nbins
                    IPSPslice(j) = sum(profiles(mei((ind==j)&pickedNeighbor,i)));
                end
                IPSPslice = IPSPslice/sum(IPSPslice)*100;
                plotyy(binranges,[IPSPslice,norm_binCounts,norm_binCounts0],binranges,binCounts*100./binCounts0);
                ylabel('%');
                legend({'Cort. Inh.','Pre. Neuron','All Neurons','Prob'},'Location','NorthWest');
                ylim([0,inf]);
                xlabel('Coeff');
                check([i;p.nv1e+bigShots],h,coeffs);
            %% theta
                subplot(ceil((nn-1)/plotPr)+1,3,9);
                binranges = 0.0:dtheta:(pi/2);
                nbins = length(binranges)-1;
                binranges(nbins+1) = pi/2+0.1;
                dThetas(~pickedNeighbor) = binranges(1)-1;
                [binCounts,ind] = histc(dThetas,binranges);
                binCounts = binCounts(1:nbins);
                mbinCounts = binCounts./sum(binCounts)*100;
                IPSPslice = zeros(nbins,1);
                mIPSPslice = zeros(nbins,1);
                for j = 1:nbins
                    %IPSPslice(j) = sum(profiles(mei((ind==j)&pickedNeighbor,i)));
                    tmp = (profiles(mei((ind==j)&pickedNeighbor,i)));
                    mIPSPslice(j) = mean(tmp);
                    IPSPslice(j) = sum(tmp);
                end
                %[hAx,h1,h2] = plotyy(dThetas(pickedNeighbor)*180/pi,EPSPs,(binranges(1:nbins)+dtheta/2)*180/pi,IPSPslice);
                IPSPslice = IPSPslice./sum(IPSPslice)*100;
                xxxx = (binranges(1:nbins)+dtheta/2)*180/pi; 
                [hAx,h1,h2] = plotyy(xxxx,mIPSPslice,xxxx,[mbinCounts,IPSPslice]);
                h1.Marker = '.';
                h2(1).LineStyle = ':';
                h2(1).Marker = 'o';
                h2(2).Marker = 's';
                ylabel(hAx(1),'IPSP (mV)');
                ylabel(hAx(2),'%');
                xlabel('\Delta\theta');
                legend({'avg IPSP','% Connection','% Total IPSP'});
                if ~isempty(format)
                    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                    print(h,['EI_RFsimilarity-',num2str(i),'-',coMatfile,'-',theme,'.',format],printDriver,dpi);
                end
            end
            mcoeff(i) = mean(coeffs(bigShots));
            scoeff(i) = std(coeffs(bigShots));
        end
        fp = 3;
        disp(['mean coeff ',num2str(mean(mcoeff),fp),'+-',num2str(std(mcoeff),fp)]);
        disp(['std of coeff ',num2str(mean(scoeff),fp),'+-',num2str(std(scoeff),fp)]);
    end
end
function cR = checkRegion(list,src,tar)
    nlist = length(list);
    cR = false(src.n,nlist);
    dx = zeros(src.n,1);
    dy = zeros(src.n,1);
    hlx = src.lx/2;
    hly = src.ly/2;
    HWHM2 =  (src.raxn + tar.rden)^2;

    r2 = tar.bound^2*HWHM2;
    for k = 1:nlist
        i = list(k);
        pick = abs(src.x-tar.x(i)) > hlx;
        dx(pick) = 2*mod(tar.x(i)+hlx,src.lx) - tar.x(i) - src.x(pick);
        %assert(sum(dx(pick) >= hlx + 1e-5)==0);
        dx(~pick) = src.x(~pick) - tar.x(i);
        %y
        pick = abs(src.y-tar.y(i)) > hly;
        dy(pick) = 2*mod(tar.y(i)+hly,src.ly) - tar.y(i) - src.y(pick);
        %assert(sum(dy(pick) >= hly + 1e-5)==0);
        dy(~pick) = src.y(~pick) - tar.y(i);

        d2 = dx.^2 + dy.^2;
        % boundary cutoff
        cR(:,k) = d2 < r2;
    end
end
function check(list,h,coeffs)
    global x y X Y LGNpos v1Map lgnStrength nSubLGN p etheta itheta rfx rfy ttt FWHM plotPr
    ymax0 = -inf;
    xmax0 = ymax0;
    ymin0 = inf;
    xmin0 = ymin0;
    figure(h);
    n = length(list);
    for k = 1:n
        i = list(k);
        if k==1
            subplot(ceil((n-1)/plotPr)+1,3,3);
        else
            ch = subplot(ceil((n-1)/plotPr)+1,plotPr,k+plotPr-1);
        end
        hold on
        subregion = length(nSubLGN{i});
        subpick = 1:subregion;
        if i <= p.nv1e
            ii = i;
            theta = etheta(i);
            sigma = p.esigma(i,subpick);
            sigmb = sigma.*p.eAspectRatio(i,subpick);
        else
            ii =i-p.nv1e;
            theta = itheta(i-p.nv1e);
            sigma = p.isigma(i-p.nv1e,subpick);
            sigmb = sigma.*p.iAspectRatio(i-p.nv1e,subpick);
        end
        if theta >=pi/2
            gtheta = theta - pi/2;
        else
            gtheta = theta + pi/2;
        end

        peak = zeros(subregion,2);
        for j = 1:subregion
            if nSubLGN{i}(j) > 0
                peak(j,:) = mean(LGNpos(v1Map{i,j},:),1)*180/pi;
            else
                peak(j,:) = mean(LGNpos(-v1Map{i,j},:),1)*180/pi;
            end
        end
        ymax = -inf;
        xmax = ymax;
        ymin = inf;
        xmin = ymin;
        xxx = zeros(length(ttt),subregion);
        yyy = xxx;
        for j = 1:subregion
            
            xxx(:,j) = rfx(peak(j,1)/180*pi,sigmb(j),sigma(j),gtheta,ttt,FWHM,FWHM)*180/pi;         
            yyy(:,j) = rfy(peak(j,2)/180*pi,sigmb(j),sigma(j),gtheta,ttt,FWHM,FWHM)*180/pi;
            ymin = min(min(min(yyy(:,j)),ymin),ymin0);
            ymax = max(max(max(yyy(:,j)),ymax),ymax0);
            xmin = min(min(min(xxx(:,j)),xmin),xmin0);
            xmax = max(max(max(xxx(:,j)),xmax),xmax0);
            if nSubLGN{i}(j) > 0
                Z = multiLGNspatialKernel(LGNpos(v1Map{i,j},:)*180/pi,nSubLGN{i}(j),lgnStrength{i,j},1);
                newZmax = max(max(Z));
                if (k==1)
                    [~,c1(j)] = contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','m','LineWidth',1);
                    [~,c2(j)] = contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','m','LineWidth',1);
                    if gtheta > pi/4
                        l(j) = plot((y-peak(j,2))*cot(gtheta)+peak(j,1),y,':m');
                    else
                        l(j) = plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':m');
                    end
                else
                    contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','r','LineWidth',1);
                    contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','r','LineWidth',1);
                    if gtheta > pi/4
                        plot((y-peak(j,2))*cot(gtheta)+peak(j,1),y,':r');
                    else
                        plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':r');
                    end
                end
            else
                Z = multiLGNspatialKernel(LGNpos(-v1Map{i,j},:)*180/pi,-nSubLGN{i}(j),lgnStrength{i,j},1);
                newZmax = max(max(Z));
                if (k==1)
                    [~,c1(j)] = contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','c','LineWidth',1);
                    [~,c2(j)] = contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','c','LineWidth',1);
                    if gtheta > pi/4
                        l(j) = plot((y-peak(j,2))*cot(gtheta)+peak(j,1),y,':c');
                    else
                        l(j) = plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':c');
                    end
                else
                    contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','b','LineWidth',1);
                    contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','b','LineWidth',1);
                    if gtheta > pi/4
                        plot((y-peak(j,2))*cot(gtheta)+peak(j,1),y,':b');
                    else
                        plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':b');
                    end
                end
            end
        end
        if k>1
            copyobj([c1,c2,l],ch);
        end
        axis equal
        axis([xmin,xmax,ymin,ymax]);
        if k>1
            %title([num2str(i),'-',num2str(subregion),' corr =',num2str(coeffs(ii),'%3.2f'), 'g\theta = ',num2str(gtheta*180/pi,'%3.2f')]);
            if i > p.nv1e
                title({[num2str(i),'-',num2str(subregion),' corr =',num2str(coeffs(ii),'%3.2f')], ['normD = ',num2str(p.inormDistance(i-p.nv1e),'%3.2f')]});
            else
                title({[num2str(i),'-',num2str(subregion),' corr =',num2str(coeffs(ii),'%3.2f')], ['normD = ',num2str(p.enormDistance(i),'%3.2f')]});
            end
        else
            title(['target neuron ',num2str(i),' g\theta = ',num2str(gtheta*180/pi,'%3.2f')]);
        end
        %axis([x(1),x(nx),y(1),y(ny)]);
        if k == 1
            ymax0 = ymax;
            xmax0 = xmax;
            ymin0 = ymin;
            xmin0 = xmin;
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
end
