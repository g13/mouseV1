function plotISI(theme,lgnfile,eps,ctheta,threads,format,ntheta)
    if nargin < 7
        ntheta = 12;
        if nargin < 6
            format = '';
        end
    end
    nbinsE = 150;
    nbinsI = 250;
    frBinsE = 15;
    frBinsI = 15;
    scale = 1000;   % 1s = 1000ms
    roundoff = 100; % ticks at 100ms
    factors = [2,5,10,20,50,100,200];

    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r150';
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
    eps = 125*2^(eps-1);
    epsStr = num2str(eps,'%04d');
    thetaStr = num2str(ctheta,'%02d');
    DIR = [theme,'/spike_wise/',epsStr,'/'];
    [tspI,l] = readSpikes(DIR,p.nv1,[thetaStr,'-spikes.dat']);
    fid = fopen([theme,'/',epsStr,'/cv.dat']);
    fread(fid,[p.nv1,4],'double');
    fr = fread(fid,[p.nv1,1],'double');
    fread(fid,[p.nv1,6],'double');
    priA = fread(fid,[p.nv1,1],'double');
    dtheta = 180/ntheta;
    ipA = round(priA*180/pi/dtheta);
    epick = [ipA(1:p.nv1e) == ctheta; false(p.nv1i,1)];
    [isiE, minE, maxE] = getISI(tspI(epick));
    ipick = [false(p.nv1e,1); ipA(p.nv1e+(1:p.nv1i)) == ctheta];
    [isiI, minI, maxI] = getISI(tspI(ipick));
    [isiMatE, edgesMatE, activeE] = binitMat(isiE,minE,maxE,nbinsE,scale,roundoff,factors);
    [isiMatI, edgesMatI, activeI] = binitMat(isiI,minI,maxI,nbinsI,scale,roundoff,factors);
    [isiVecE, edgesVecE] = binitVec(isiE,minE,maxE,nbinsE,scale,roundoff,factors);
    [isiVecI, edgesVecI] = binitVec(isiI,minI,maxI,nbinsI,scale,roundoff,factors);
    h = figure;
    ax = subplot(2,2,1);
    plotsubMat(isiMatE,edgesMatE,scale,ax);
    xlabel('ISI (ms)');
    ylabel('pdf %');
    ylim([0,inf]);
    ax = subplot(2,2,2);
    plotsubMat(isiMatI,edgesMatI,scale,ax);
    xlabel('ISI (ms)');
    ylabel('pdf %');
    ylim([0,inf]);
    ax = subplot(2,2,3);
    plotsubVec(isiVecE,edgesVecE,scale,ax);
    xlabel('ISI (ms)');
    ylabel('pdf %');
    ylim([0,inf]);
    ax = subplot(2,2,4);
    plotsubVec(isiVecI,edgesVecI,scale,ax);
    xlabel('ISI (ms)');
    ylabel('pdf %');
    ylim([0,inf]);
    if ~isempty(format)
        fname = [theme,'/ISI-',epsStr,'-',thetaStr,'-',theme,'.',format];
        if strcmp(format,'fig')
            saveas(h,fname);
        else
            print(h,fname,printDriver,dpi);
        end
    end
    h = figure;
    %subplot(2,2,1);
    efr = fr(epick);
    [H,frEdgesE,binIDE] = histcounts(efr(activeE),frBinsE);
    %edges = frEdgesE(1:frBinsE)- (frEdgesE(2)-frEdgesE(1))/2;
    %bar(edges,H);
    ax = subplot(2,1,1);
    heatFrISI(isiMatE,frEdgesE,binIDE,edgesMatE,scale,ax);
    xlabel('ISI (ms)');
    ylabel('FR (Hz)');
    title('Excitatory');
    %subplot(2,2,3);
    ifr = fr(ipick);
    [H,frEdgesI,binIDI] = histcounts(ifr(activeI),frBinsI);
    %edges = frEdgesI(1:frBinsI)- (frEdgesI(2)-frEdgesI(1))/2;
    %bar(edges,H);
    ax = subplot(2,1,2);
    heatFrISI(isiMatI,frEdgesI,binIDI,edgesMatI,scale,ax);
    xlabel('ISI (ms)');
    ylabel('FR (Hz)');
    title('Inhibitory');
    if ~isempty(format)
        fname = [theme,'/ISIheat-',epsStr,'-',thetaStr,'-',theme,'.',format];
        if strcmp(format,'fig')
            saveas(h,fname);
        else
            print(h,fname,printDriver,dpi);
        end
    end
end
function [isi, minISI, maxISI] = getISI(tspI)
    n = length(tspI);
    isi = cell(1,n);
    maxtmp = zeros(n,1);
    mintmp = zeros(n,1)+inf;
    parfor i=1:n
    %for i = 1:n
        isi{i} = diff(tspI(i).tsp);
        assert(sum(isi{i}>0) == length(isi{i}));
        if ~isempty(isi{i})
            maxtmp(i) = max(isi{i});
            mintmp(i) = min(isi{i});
        end
    end
    maxISI = max(maxtmp);
    minISI = min(mintmp);
end
function [isiMat, edges, active] = binitMat(isi,minISI,maxISI,nbins,scale,roundoff,factors)
    n = length(isi);
    minISI = floor(floor(minISI*scale)/roundoff)*roundoff;
    maxISI = ceil(ceil(maxISI*scale)/roundoff)*roundoff;
    len = maxISI-minISI
    dISI = len/nbins;
    [~,facID] = min(abs(dISI - factors));
    dISI = factors(facID);
    edges = (minISI:dISI:maxISI)/scale;
    nbins = length(edges)-1
    isiMat = zeros(n,nbins);
    active = true(n,1);
    parfor i=1:n
    %for i = 1:n
        if ~isempty(isi{i})
            [tmp, ~] = histcounts(isi{i},edges);
            isiMat(i,:) = tmp;
        else
            active(i) = false;
        end
    end
    isiMat = isiMat(active,:);
end
function [isiVec, edges] = binitVec(isi,minISI,maxISI,nbins,scale,roundoff,factors)
    n = length(isi);
    minISI = floor(floor(minISI*scale)/roundoff)*roundoff;
    maxISI = ceil(ceil(maxISI*scale)/roundoff)*roundoff;
    dISI = ((maxISI-minISI)/nbins);
    [~,facID] = min(abs(diff(dISI - (roundoff./factors))));
    dISI = roundoff/factors(facID); 
    edges = (minISI:dISI:maxISI)/scale;
    nbins = length(edges)-1;
    [isiVec, ~] = histcounts(cell2mat(isi),edges);
end
function plotsubMat(isiMat,edges,scale,ax)
    isiMat_pdf = isiMat*100./(sum(isiMat,2)*ones(1,size(isiMat,2)));
    m = mean(isiMat_pdf);
    l = m-quantile(isiMat_pdf,0.25);
    u = quantile(isiMat_pdf,0.75) - m;
    edges = edges(1:length(edges)-1) + (edges(2) - edges(1))/2;
    axes(ax);
    errorbar(ax,edges*scale, m, l, u);
end
function plotsubVec(isiVec,edges,scale,ax)
    s = sum(isiVec);
    if s > 0
        isiVec = isiVec/s;
    end
    edges = edges(1:length(edges)-1) + (edges(2) - edges(1))/2;
    axes(ax);
    bar(ax,edges*scale,isiVec);
end
function heatFrISI(isiMat,frEdges,binID,edges,scale,ax)
    nFrBins = length(frEdges)-1;
    nbins = length(edges)-1;
    dataPair = zeros(nbins,nFrBins);
    axes(ax);
    for i=1:nFrBins
        if sum(binID==i) > 0
            dataPair(:,i) = sum(isiMat(binID==i,:),1)';
        end
    end
    imagesc([1,nbins],[nFrBins,1],dataPair');
    %daspect([1,1,1]);
    colormap('hot');
    %axis tight
    %box on

    dedges = (edges(2) - edges(1)) * scale;
    xtick = get(gca,'XTick');
    xticklabel = num2str(edges(1) + (xtick'-1)*dedges);
    xtick = xtick-0.5;
    set(gca,'XTick',xtick,'XTickLabel',xticklabel);

    dFrEdges = (frEdges(2) - frEdges(1));
    ytick = get(gca,'YTick');
    yticklabel = num2str(frEdges(1) + (ytick'-1)*dFrEdges);
    ytick = ytick-0.5;
    set(gca,'YTick',ytick,'YTickLabel',flipud(yticklabel));
end
