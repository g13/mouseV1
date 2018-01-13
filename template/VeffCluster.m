function VeffCluster(theme,eps,angle,clusterFile,twindow,t1,format,mdt,dt)
    ntotal = 10800;
    if nargin < 9
        dt = 1e-3;
        if nargin < 8
            mdt = 5;
            if nargin < 7
                format = '';
            end
        end 
    end
    t0 = t1-twindow;  %s
    t = t0:dt:t1;
    % left bot wid height
    FontSize = 12;
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
    FontSize = 10;
    epsfdr = num2str(eps,'%04d');
    anglefdr = num2str(angle,'%02d');

    fid = fopen([theme,'/sample_id.dat']);
    sampleID = fread(fid,[1,inf], 'int');
    nsample = length(sampleID);
    disp(['nsample = ',num2str(nsample)]);
    fclose(fid);

    filepath = [theme,'/samples/',epsfdr,'/',anglefdr,'-samples.dat'];
    disp(filepath);
    fid = fopen(filepath);
    data = fread(fid,[8*nsample,inf],'double');
    ndt = size(data,2);
    assert(ndt==length(t));
    data = reshape(data,[8,nsample,ndt]);
    fclose(fid);
    load(clusterFile);
    Veff = reshape(data(7,:,:),[nsample,ndt]);
    Veff = Veff';
    ne = length(clusterID);
    Veff = Veff(:,1:ne);
    if mdt > 1
        rmdt = mod(ndt,mdt);
        ndiv = floor((ndt-rmdt)/mdt);
        if rmdt > 0
            pVeff = reshape(Veff(1:ndt-rmdt,:),[mdt,ndiv,ne]);
            pVeff = mean(pVeff);
            pVeff = reshape(pVeff,[ndiv,ne]);
            pVeff = [pVeff;mean(Veff(ndt-rmdt+1:ndt,:),1)];
            cVeff = zeros(ndiv+1,nCluster);
        else
            pVeff = reshape(Veff,[mdt,ndiv,ne]);
            pVeff = mean(pVeff);
            pVeff = reshape(pVeff,[ndiv,ne]);
            cVeff = zeros(ndiv,nCluster);
        end
    else
        pVeff = Veff;
        cVeff = zeros(ndiv,nCluster);
    end
    clear Veff;
    if ntotal == nsample
        [~,clusterSeq] = sort(clusterID,'ascend');
        pVeff = pVeff(:,clusterSeq);
        j = 0;
        for i = 1:nCluster
            cVeff(:,i) = mean(pVeff(:,j+(1:clusterSize(i))),2); 
            j = j + clusterSize(i);
        end
        cVeff = cVeff./max(cVeff(:));
        sCluster = nCluster;
    else
        pickCluster = zeros(nCluster,1);
        for ii=1:ne
            i = sampleID(ii);
            cVeff(:,clusterID(i)) = cVeff(:,clusterID(i)) + pVeff(:,i);
            pickCluster(clusterID(i)) = pickCluster(clusterID(i))+1;
        end
        picked = pickCluster>0;
        cVeff = cVeff(:,picked)./pickCluster(picked);
        cVeff = cVeff./max(cVeff(:));
        sCluster = sum(picked);
    end
    h = figure;
    [tick,n0] = autoAxis(t0,t1,11,[t0,t1]);
    tickPosX = linspace(0.5,ndiv-1+0.5,n0);
    tickLabelX = num2str(tick');

    [tick,n0] = autoAxis(0,sCluster,min(11,sCluster),[0,sCluster]);
    tickPosY = linspace(0.5,sCluster-1+0.5,n0);
    tickLabelY = num2str(tick');

    colormap(redOnly);
    imagesc([1,ndiv],[sCluster,1],cVeff);
    set(gca,'YTickLabel',flipud(tickLabelY),'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
    xlabel('time /s');
    ylabel('Cluster ID');
    colorbar;
    title(['averaged over ',num2str(mdt),'ms and cluster']);
    
    if ~isempty(format)
        figname = ['cVeff-',theme,'-',epsfdr,'-',anglefdr,'_',num2str(nthres),'-',num2str(cthres),'.',format];
        if strcmp(format,'fig')
            savefig(h,figname);
        else
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(h,figname,printDriver,dpi);
        end
    end

end
