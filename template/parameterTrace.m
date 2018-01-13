function parameterTrace(fdrlist,pLabel,pTickLabel,lgnfile,format,theme)
    pPosition = [0, 0, 1280, 720];
    set(0, 'DefaultFigurePosition', pPosition);
    LineWidth = 2;
    set(groot,'defaultLineLineWidth',LineWidth);
    set(groot,'defaultErrorbarLineWidth',LineWidth);
    FontSize = 14;
    set(groot,'defaultAxesFontSize',FontSize);
    set(groot,'defaultTextFontSize',FontSize);
    LegendOffset = 1;
    set(groot,'defaultLegendFontSize',FontSize-LegendOffset);
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


    nfdr = length(fdrlist);
    load(lgnfile);
    cvE     = zeros(nfdr,4);
    cvI     = zeros(nfdr,4);
    %gcv    = zeros(nfdr,2);
    scE     = zeros(nfdr,3);
    scI     = zeros(nfdr,3);
    spontE  = zeros(nfdr,2);
    evokedE = zeros(nfdr,2);
    spontI  = zeros(nfdr,2);
    evokedI = zeros(nfdr,2);
    widthE  = zeros(nfdr,2);
    condR = zeros(nfdr,4);
    ep = 1:p.nv1e;
    ip = p.nv1e + (1:p.nv1i);
    for i = 1:nfdr
        fdr = fdrlist{i};
        load([fdr,'-tcData-x4.mat'],'nP','ipC','ioC')

        target = 1-nP(2).cv(ep);
        cvE(i,1) = mean(target);
        cvE(i,2) = std(target);
        target = 1-nP(4).cv(ep);
        cvE(i,3) = mean(target);
        cvE(i,4) = std(target);

        target = 1-nP(2).cv(ip);
        cvI(i,1) = mean(target);
        cvI(i,2) = std(target);
        target = 1-nP(4).cv(ip);
        cvI(i,3) = mean(target);
        cvI(i,4) = std(target);
        
        target = nP(4).sc(ep);
        edges = 0.0:0.2:2.0;
        [N,~] = histcounts(target,edges);
        [~,ind] = max(N);
        scE(i,1) = (ind-1)*0.2+0.1;
        scE(i,2) = mean(target);
        scE(i,3) = std(target);

        target = nP(4).sc(ip);
        [N,~] = histcounts(target,edges);
        [~,ind] = max(N);
        scI(i,1) = (ind-1)*0.2+0.1;
        scI(i,2) = mean(target);
        scI(i,3) = std(target);

        target = nP(1).br(ep);
        spontE(i,1) = mean(target);
        spontE(i,2) = std(target);

        target = nP(1).br(ip);
        spontI(i,1) = mean(target);
        spontI(i,2) = std(target);

        target = nP(4).pkrate(ep);
        evokedE(i,1) = mean(target);
        evokedE(i,2) = std(target);

        target = nP(4).pkrate(ip);
        evokedI(i,1) = mean(target);
        evokedI(i,2) = std(target);

        target = mean(ipC(2).gI(ep,:),2) ./ mean(ipC(2).gE(ep,:)+ipC(2).gLGN(ep,:),2);
        condR(i,1) = mean(target);
        target = mean(ioC(2).gI(ep,:),2) ./ mean(ioC(2).gE(ep,:)+ioC(2).gLGN(ep,:),2);
        condR(i,2) = mean(target);

        target = mean(ipC(4).gI(ep,:),2) ./ mean(ipC(4).gE(ep,:)+ipC(4).gLGN(ep,:),2);
        condR(i,3) = mean(target);
        target = mean(ioC(4).gI(ep,:),2) ./ mean(ioC(4).gE(ep,:)+ioC(4).gLGN(ep,:),2);
        condR(i,4) = mean(target);

        clear nP ioC ipC
        fitfile = [theme,'/',theme,'-4fitted.mat'];
        if exist(fitfile,'file');
            load(fitfile);
            target = sigfsq2(ep) * sqrt(log(2));
            width(i,1) = mean(target);
            width(i,2) = std(target);
        end
    end
    hCV = figure;
    subplot(1,2,1);
    hold on
    x = 0:1;
    y = 0:1;
    plot(x,y,':k');
    for i=1:nfdr
        %if nfdr > 1
        %    col = [0,0.3 + 0.7*(i-1)/(nfdr-1),1];
        %else
        %    col = [0,0.3,1];
        %end
        if nfdr > 1
            col = [0,1,0.3 + 0.7*(i-1)/(nfdr-1)];
        else
            col = [0,1,1];
        end
        col = hsv2rgb(col);
        h(i) = errorbar(cvE(i,1),cvE(i,3),cvE(i,4),'.','MarkerSize',32,'Color',col);
        x = linspace(cvE(i,1)-cvE(i,2),cvE(i,1)+cvE(i,2),2);
        y = cvE(i,3)*ones(length(x),1);
        plot(x,y,'-','Color',col);
    end
    legend(h,pTickLabel,'Location','SouthEast');
    xlabel('gOSI (25%)');
    ylabel('gOSI (100%)');
    title('Exc');

    subplot(1,2,2);
    hold on
    x = 0:1;
    y = 0:1;
    plot(x,y,':k');
    for i=1:nfdr
        %if nfdr > 1
        %    col = [2/3,0.3+0.7*(i-1)/(nfdr-1),1];
        %else
        %    col = [2/3,0.3,1];
        %end
        if nfdr > 1
            col = [2/3,1,0.3 + 0.7*(i-1)/(nfdr-1)];
        else
            col = [2/3,1,1];
        end
        col = hsv2rgb(col);
        h(i) = errorbar(cvI(i,1),cvI(i,3),cvI(i,4),'.','MarkerSize',32,'Color',col);
        x = linspace(cvI(i,1)-cvI(i,2),cvI(i,1)+cvI(i,2),2);
        y = cvI(i,3)*ones(length(x),1);
        plot(x,y,'-','Color',col);
    end
    legend(h,pTickLabel,'Location','NorthWest');
    xlabel('gOSI (25%)');
    ylabel('gOSI (100%)');
    title('Inh')
    ylim([0,0.5]);
    xlim([0,0.5]);

    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hCV,['CV-',theme,'.',format]);
        else
            print(hCV,['CV-',theme,'.',format],printDriver,dpi);
        end
    end
    hSC = figure;
    subplot(1,2,1);
    hold on
    errorbar(1:nfdr,scE(:,2),scE(:,3),'.r','MarkerSize',32);
    plot(1:nfdr,scE(:,1),'*r');
    xlabel(pLabel);
    set(gca,'XTickLabel',pTickLabel,'XTick',1:nfdr);
    ylabel('F_1/F_0')
    ylim([0,2]);
    title('Exc');
    
    subplot(1,2,2);
    hold on
    errorbar(1:nfdr,scI(:,2),scI(:,3),'.b','MarkerSize',32);
    plot(1:nfdr,scI(:,1),'*b');
    xlabel(pLabel);
    set(gca,'XTickLabel',pTickLabel,'XTick',1:nfdr);
    ylabel('F_1/F_0')
    ylim([0,2]);
    title('Inh');
    
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hSC,['SC-',theme,'.',format]);
        else
            print(hSC,['SC-',theme,'.',format],printDriver,dpi);
        end
    end
    hFR = figure;
    subplot(1,2,1);
    hold on
    target = [zeros(nfdr,1),spontE(:,1)];
    b = bar([0,1],target');
    for i = 1:nfdr
        if nfdr > 1
            b(i).FaceColor = hsv2rgb([0,1.0,0.3+0.7*(i-1)/(nfdr-1)]);
        else
            b(i).FaceColor = hsv2rgb([0,1.0,1.0]);
        end
    end
    legend(pTickLabel,'Location','NorthWest');
    
    target = [spontI(:,1),zeros(nfdr,1)];
    b = bar([2,3],target');
    for i = 1:nfdr
        if nfdr > 1
            b(i).FaceColor = hsv2rgb([2/3,1.0,0.3+0.7*(i-1)/(nfdr-1)]);
        else
            b(i).FaceColor = hsv2rgb([2/3,1.0,1.0]);
        end
    end

    l = 0.68;
    shift = 0.16;
    dw = l/(nfdr+1);
    xe(1:nfdr) = 0.5 + shift + linspace(dw/2,l-dw/2,nfdr);
    xi(1:nfdr) = 1.5 + shift + linspace(dw/2,l-dw/2,nfdr);
    targetstd = spontE(:,2)';
    target = spontE(:,1)';
    errorbar(xe,target,zeros(length(target),1),targetstd,'.k');

    targetstd = spontI(:,2)';
    target = spontI(:,1)';
    errorbar(xi,target,zeros(length(target),1),targetstd,'.k');

    xlim([0.5,2.5]);
    set(gca,'XTickLabel',{'Exc','Inh'},'XTick',1:2);
    ylabel('Spont. Rate /Hz');

    subplot(1,2,2);
    hold on
    target = [zeros(nfdr,1),evokedE(:,1)];
    b = bar([0,1],target');
    for i = 1:nfdr
        if nfdr > 1
            b(i).FaceColor = hsv2rgb([0,1.0,0.3+0.7*(i-1)/(nfdr-1)]);
        else
            b(i).FaceColor = hsv2rgb([0,1.0,1.0]);
        end
    end
    legend(pTickLabel,'Location','NorthWest');

    target = [evokedI(:,1),zeros(nfdr,1)];
    b = bar([2,3],target');
    for i = 1:nfdr
        if nfdr > 1
            b(i).FaceColor = hsv2rgb([2/3,1.0,0.3+0.7*(i-1)/(nfdr-1)]);
        else
            b(i).FaceColor = hsv2rgb([2/3,1.0,1.0]);
        end
    end

    dw = l/(nfdr+1);
    xe = 0.5 + shift + linspace(dw/2,l-dw/2,nfdr);
    xi = 1.5 + shift + linspace(dw/2,l-dw/2,nfdr);
    targetstd = evokedE(:,2)';
    target = evokedE(:,1)';
    errorbar(xe,target,zeros(length(target),1),targetstd,'.k');

    targetstd = evokedI(:,2)';
    target = evokedI(:,1)';
    errorbar(xi,target,zeros(length(target),1),targetstd,'.k');

    xlim([0.5,2.5]);
    set(gca,'XTickLabel',{'Exc','Inh'},'XTick',1:2);
    ylabel('Evoked Rate /Hz');

    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hFR,['FR-',theme,'.',format]);
        else
            print(hFR,['FR-',theme,'.',format],printDriver,dpi);
        end
    end
    hcondR = figure;
    yy = zeros(nfdr,1);
    for i=1:nfdr
        if nfdr > 1
            col1 = [0,0.35,0.3+0.7*(i-1)/(nfdr-1)];
        else
            col1 = [0,0.35,1.0];
        end
        col2 = col1;
        col2(2) = 1;
        col1 = hsv2rgb(col1);
        col2 = hsv2rgb(col2);
        ax(i)= subplot(2,nfdr,i);
        hold on
        plot(1:2,[condR(i,2),condR(i,1)],'-*','Color',col1);
        plot(1:2,[condR(i,4),condR(i,3)],'-*','Color',col2);
        legend({'25%','100%'});
        set(gca,'XTicklabel',{'orth.','pref.'},'XTick',1:2);
        xlim([0.5,2.5]);
        ytmp = ylim;
        yy(i) = ytmp(2);
        ylabel('g_I:g_E');
        title(pTickLabel(i));
    end
    [~,id] = sort(yy,'descend');
    ylim(ax(id(1)),[0,yy(i)*1.3]);
    linkaxes(ax(id),'xy');
    subplot(2,nfdr,nfdr+1);
    hold on
    for i=1:nfdr
        %if nfdr > 1
        %    col1 = [0+1/6*(i-1)/(nfdr-1),0.40,0.9];
        %else
        %    col1 = [0,0.40,0.9];
        %end
        if nfdr > 1
            col1 = [0,0.35,0.3+0.7*(i-1)/(nfdr-1)];
        else
            col1 = [0,0.35,1.0];
        end
        col2 = col1;
        col2(2) = 1;
        col1 = hsv2rgb(col1);
        col2 = hsv2rgb(col2);
        plot(1:2,[condR(i,2),condR(i,1)],'-*','Color',col1);
        h(i) = plot(1:2,[condR(i,4),condR(i,3)],'-*','Color',col2);
        %if i == 1
        %    legend({'25%','100%'},'Location','NorthWest');
        %end
    end
    set(gca,'XTicklabel',{'orth.','pref.'},'XTick',1:2);
    ylabel('g_I:g_E');
    xlim([0.5,2.5]);
    yy = ylim;
    ylim([0,yy(2)*1.3]);
    ah=axes('position',get(gca,'position'),...
                'visible','off');
    legend(ah,h,pTickLabel,'Location','NorthEast');
%    subplot(2,1,2);
%    target = [(condR(:,4)./condR(:,3))./(condR(:,2)./condR(:,1))];
%    bar(target');
%    title('100%(pref/orth) : 25%(pref/orth) of g_I:g_E');
%    set(gca,'XTicklabel',pTickLabel,'XTick',1:nfdr);

    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hcondR,['condR-',theme,'.',format]);
        else
            print(hcondR,['condR-',theme,'.',format],printDriver,dpi);
        end
    end
    hWidth = figure;
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(format,'fig')
            saveas(hWidth,['Width-',theme,'.',format]);
        else
            print(hWidth,['Width-',theme,'.',format],printDriver,dpi);
        end
    end
end
