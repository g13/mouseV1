function plotsample(theme,eps,angle,nrange,sampleID,nsample,t,nlgn,period,format,processed,lgnfile,tshift,nd)
    ntotal=10800;
    FontSize = 10;
    vTheta = 4.375;
    pPosition = [0, 0, 1280, 720];
    set(0,'DefaultAxesFontSize',FontSize)
    set(0, 'DefaultFigurePosition', pPosition);
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r300';
    end
    % left bot wid height
    fid = fopen([theme,'/','cMatrix_summary']);
    discard = fread(fid, [2,ntotal],'int32');
    preE = fread(fid, [1,ntotal],'int32');
    preI = fread(fid, [1,ntotal],'int32');
    clear discard
    mat = fread(fid, [ntotal,ntotal],'int8');
    fclose(fid);

    epsfdr = num2str(eps*1000,'%04d');
    anglefdr = num2str(angle,'%02d');
    if processed
        filepath = [theme,'/samples/',epsfdr,'/',anglefdr,'-samples.dat'];
    else
        filepath = [theme,'/',epsfdr,'/',anglefdr,'/samples.dat'];
    end
    disp(filepath);
    fid = fopen(filepath);
    data = fread(fid,[9*nsample,inf],'double');
    fclose(fid);
    ndt = size(data,2);
    disp(['ndt = ',num2str(ndt)]);
    data = reshape(data,[9,nsample,ndt]);
    trange = 1:ndt;
    periods = t(trange(1)):period:t(trange(end));
    marks = zeros(1,length(periods));

    rasterfile = [theme,'-',epsfdr,'-',anglefdr,'-rasterStruct.mat'];
    if exist(rasterfile,'file')
        load(rasterfile);
    else
        rasterSave(theme,lgnfile,eps*1000,angle,processed)
        load(rasterfile);
    end 
    for i=1:length(nrange)
        isample = reshape(data(:,nrange(i),:),[9,ndt]);
        isample = isample';
        h = figure;
            subplot(6,1,1);hold on;grid on;
        plot(t(trange),isample(trange,1));
        plot(periods,marks,'*r');
        plot(periods,marks,':g');
        ylabel('v');
        xlim([t(trange(1)),t(trange(end))]);
        ylim([-inf,vTheta]);
        if angle~= 12
        title({['neuron ',num2str(sampleID(i)), ' nLGN = ',num2str(nlgn(i)),...
            ' eps =', num2str(eps*100),'%',...
            ' \theta = ',num2str(angle/16*180),'^{o}'],['nD =', num2str(nd(i))]});
        else
            title({['neuron ',num2str(sampleID(i)), ' nLGN = ',num2str(nlgn(i)),...
            ' eps =', num2str(eps*100),'%',...
            ' spontaneous'],['nD =', num2str(nd(i))]});
        end
            subplot(6,1,2);hold on;grid on;
        plot(t(trange),isample(trange,2));
        plot(periods,marks,'*r');
        ylabel('gLGN');
        xlim([t(trange(1)),t(trange(end))]);
            subplot(6,1,3);hold on;grid on;
        plot(t(trange),isample(trange,3));
        plot(t(trange),isample(trange,5),'m');
        plot(periods,marks,'*r');
        iii = sampleID(i);
        if iii <= 8640
            maxge = max(isample(trange,3));
            [~, id] = sort(mat(iii,:),'descend');
            id = id(1:preE(iii));
            rsdata = ones(sum(l(id)),1);
            last = 0;
            for j=1:preE(iii)
               rsdata(last+(1:l(id(j)))) = tspI(id(j)).tsp;
               last = last + l(id(j));
            end
            N = histcounts(rsdata,t(trange)+tshift);
            plot(t(trange(1:end-1)),N./max(N)*maxge,'g');
            np = 10;
            for ii=1:np
                pick = tspI(id(ii)).tsp < (t(trange(end))+tshift) & tspI(id(ii)).tsp > (t(trange(1))+tshift);
                lp = sum(pick);
                dot = (np-(ii-1))/np * maxge * ones(lp,1);
                plot(tspI(id(ii)).tsp(pick)-tshift,dot,'.');
            end
        end
        ylabel('ge');
        xlim([t(trange(1)),t(trange(end))]);

            subplot(6,1,4);hold on;grid on;
        [ax,~,~] = plotyy(t(trange),isample(trange,4),t(trange),isample(trange,9));
        plot(t(trange),isample(trange,6),'m');
        plot(periods,marks,'*r');
        ylabel('gi & gp');
        xlim([t(trange(1)),t(trange(end))]);
        ylim(ax(2),[0,inf]);

            subplot(6,1,5);hold on;grid on;
        plot(t(trange),isample(trange,7));
        plot(periods,marks,'*r');
        ylabel('gtot');
        xlim([t(trange(1)),t(trange(end))]);
            subplot(6,1,6);hold on;grid on;
        plot(t(trange),isample(trange,8)./isample(trange,7));
        plot(periods,marks,'*r');
        ylabel('vs');
        axis tight
        xlabel('t /10ms');
        if ~isempty(format)
            figname = [theme,'/sample-',num2str(sampleID(i)),'-',theme,'-',epsfdr,'-',anglefdr,'.',format];
            if strcmp(format,'fig');
                savefig(h,figname);
            else
                set(gcf,'Renderer','Painters')
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
                print(h,figname,printDriver,dpi);
            end
            disp([figname,' written']);
        end
    end
end
