function CondVsAxis(pfile,theme,ntheta,loadData,format,dir)
    if nargin < 5
        dir = theme;
        if nargin < 4
            format = '';
            if nargin < 3
                loadData = false;
                if nargin < 2
                    ntheta = 12;
                end
            end
        end
    end
    load(pfile,'p');
    set(gca,'YColor','none');
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r300';
    end
    dir = [dir,'/'];
    contrastRange = [2,4];
    vE = 2.8;
    vI = -0.4;
    vT = 1.0;
    vL = 0;
    lw = 2;
    if ~loadData || ~exist([theme,'-Cond.mat']);
        gEvE = zeros(p.nv1,2,length(contrastRange));
        gEf0f1vE = gEvE;
        gIvI = gEvE;
        gIf0f1vI = gEvE;
        gLGNvE = gEvE;
        gLGNf0f1vE = gEvE;
        Vs = gEvE;
        Vsf0f1 = gEvE;
        for ii = 1:length(contrastRange)
            i = contrastRange(ii);
            disp(['loading c',num2str(i)]);
            eps = 125*2^(i-1);
            conLabel(i) = {[num2str(eps/10,'%.1f'),'%']};
            DIR = [theme,'/',num2str(eps,'%04d')];
            [ipA,ioA,gtot,ge,gesc,gi,glgn,glgnsc,vs,vssc] = readCond(DIR,ntheta,p.nv1);
            disp('loaded, recalibrating');
            for k=1:p.nv1
                gT = gtot(k,ipA(k));
                gTf0f1 = gtot(k,ipA(k)) + ge(k,ipA(k))*gesc(k,ipA(k)) + glgn(k,ipA(k))*glgnsc(k,ipA(k));

                gEvE(k,1,ii) = ge(k,ipA(k))*vE/gT;
                gEf0f1vE(k,1,ii) = ge(k,ipA(k)).*(1+gesc(k,ipA(k)))*vE/gTf0f1;

                gLGNvE(k,1,ii) = glgn(k,ipA(k))*vE/gT;
                gLGNf0f1vE(k,1,ii) = glgn(k,ipA(k)).*(1+glgnsc(k,ipA(k)))*vE/gTf0f1;

                gIvI(k,1,ii) = gi(k,ipA(k))*vI/gT;
                gIf0f1vI(k,1,ii) = gi(k,ipA(k))*vI/gTf0f1;

                Vs(k,1,ii) = vs(k,ipA(k));
                Vsf0f1(k,1,ii) = vs(k,ipA(k)).*(1+vssc(k,ipA(k)));

                gT = gtot(k,ioA(k));
                gTf0f1 = gtot(k,ioA(k)) + ge(k,ioA(k))*gesc(k,ioA(k)) + glgn(k,ioA(k))*glgnsc(k,ioA(k));

                gEvE(k,2,ii) = ge(k,ioA(k))*vE/gT;
                gEf0f1vE(k,2,ii) = ge(k,ioA(k)).*(1+gesc(k,ioA(k)))*vE/gTf0f1;
                gLGNvE(k,2,ii) = glgn(k,ioA(k))*vE/gT;
                gLGNf0f1vE(k,2,ii) = glgn(k,ioA(k)).*(1+glgnsc(k,ioA(k)))*vE/gTf0f1;
                gIvI(k,2,ii) = gi(k,ioA(k))*vI/gT;
                gIf0f1vI(k,2,ii) = gi(k,ioA(k))*vI/gTf0f1;

                Vs(k,2,ii) = vs(k,ioA(k));
                Vsf0f1(k,2,ii) = vs(k,ioA(k)).*(1+vssc(k,ioA(k)));
            end
            disp('recalibrated');
        end
        save([theme,'-Cond.mat'],'gLGNvE','gIvI','gEvE','Vs','gLGNf0f1vE','gEf0f1vE','gIf0f1vI','Vsf0f1','-v7.3');
    else
        load([theme,'-Cond.mat']);
    end
    epick = 1:p.nv1e;
    ipick = p.nv1e+(1:p.nv1i);
    pick = epick;
    h = figure;
    s = [0.6,1.0];
    LineSpec = {'-',':'};
    for ii = 1:length(contrastRange)
        ax = subplot(2,1,ii);
        hold on
        plot(vE,0,'.k','MarkerSize',12);
        plot(vI,0,'.k','MarkerSize',12);
        plot(vL,0,'.k','MarkerSize',12);
        yl = [-0.05,0.05];
        c = [0,s(ii),1];
        for p = 1:2
            plot(ones(1,2)*mean(gEvE(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        c = [2/3,s(ii),1];
        for p = 1:2
            plot(ones(1,2)*mean(gIvI(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        c = [1/3,s(ii),1];
        for p = 1:2
            plot(ones(1,2)*mean(gLGNvE(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        c = [0,0,1-s(ii)];
        for p = 1:2
            plot(ones(1,2)*mean(Vs(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        pseudoVs = mean(gIvI(pick,:,ii))+mean(gEvE(pick,:,ii))+mean(gLGNvE(pick,:,ii));
        c = [0,0,1-s(ii)];
        for p = 1:2
            plot(ones(1,2)*pseudoVs(p),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw-1);
        end
        
        ylim(yl);
        xlim([-0.5,1.3]);
        title([num2str(osi(mean(Vs(pick,1,ii)),mean(Vs(pick,2,ii)))),',',num2str(osi(pseudoVs(1),pseudoVs(2)))]);
        box off
        set(ax,'YColor','none','Color','none','YTick',[],'YTickLabel','');
        ax.XAxisLocation = 'origin';
        axis equal
        set(ax,'XTick',[vI,vL,vT],'XTickLabel',{'V_{I}','V_{L}','V_{T}'});
    end
    if ~isempty(format)
        fname = [theme,'/CondAxisF0-',theme];
        saveas(h,fname);
        print(h,[fname,'.',format],printDriver,dpi);
    end
    h = figure;
    for ii = 1:length(contrastRange)
        ax = subplot(2,1,ii);
        hold on
        plot(vE,0,'.k','MarkerSize',12);
        plot(vI,0,'.k','MarkerSize',12);
        plot(vL,0,'.k','MarkerSize',12);
        yl = [-0.05,0.05];
        c = [0,s(ii),1];
        for p = 1:2
            plot(ones(1,2)*mean(gEf0f1vE(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        c = [2/3,s(ii),1];
        for p = 1:2
            plot(ones(1,2)*mean(gIvI(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        c = [1/3,s(ii),1];
        for p = 1:2
            plot(ones(1,2)*mean(gLGNf0f1vE(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        c = [0,0,1-s(ii)];
        for p = 1:2
            plot(ones(1,2)*mean(Vsf0f1(pick,p,ii)),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw);
        end
        pseudoVsf0f1 = mean(gIf0f1vI(pick,:,ii))+mean(gEf0f1vE(pick,:,ii))+mean(gLGNf0f1vE(pick,:,ii));
        c = [0,0,1-s(ii)];
        for p = 1:2
            plot(ones(1,2)*pseudoVsf0f1(p),yl,LineSpec{p},'Color',hsv2rgb(c),'LineWidth',lw-1);
        end
        
        ylim(yl);
        xlim([-0.5,1.3]);
        box off
        set(ax,'YColor','none','Color','none','YTick',[],'YTickLabel','');
        set(ax,'XTick',[vI,vL,vT],'XTickLabel',{'V_{I}','V_{L}','V_{T}'});
        ax.XAxisLocation = 'origin';
        title([num2str(osi(mean(Vsf0f1(pick,1,ii)),mean(Vsf0f1(pick,2,ii)))),',',num2str(osi(pseudoVsf0f1(1),pseudoVsf0f1(2)))]);
        axis equal
    end
    if ~isempty(format)
        fname = [theme,'/CondAxisF0F1-',theme];
        saveas(h,fname);
        print(h,[fname,'.',format],printDriver,dpi);
    end
end
