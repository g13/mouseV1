function plotISI(theme,lgnfile,eps,ctheta,threads,format)
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r100';
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
    disp(size(tspI));
    disp(p.nv1e);
    tmp = tspI(1:p.nv1e);
    [isiE, minE, maxE] = getISI(tmp);
    [isiI, minI, maxI] = getISI(tspI(p.nv1e+(1:p.nv1i)));
    nbins = 15;
    scale = 1000;
    [isiMatE, edgesE] = binit(isiE,minE,maxE,nbins,scale);
    [isiMatI, edgesI] = binit(isiI,minI,maxI,nbins,scale);
    h = figure;
    subplot(1,2,1)
    plotsub(isiMatE,edgesE,scale);
    subplot(1,2,2)
    plotsub(isiMatI,edgesI,scale);
    xlabel('ISI (ms)');
    ylabel('pdf %');
    if ~isempty(format)
        fname = [theme,'/ISI-',epsStr,'-',thetaStr,'-',theme,'.',format];
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
        isi{i} = diff(tspI(i).tsp);
        maxtmp(i) = max(isi{i});
        mintmp(i) = min(isi{i});
    end
    maxISI = max(maxtmp);
    minISI = min(mintmp);
end
function [isiMat, edges] = binit(isi,minISI,maxISI,nbins,scale)
    n = length(isi);
    isiMat = zeros(n,nbins);
    minISI = floor(minISI*scale);
    maxISI = ceil(maxISI*scale);
    lhist = maxISI-minISI;
    edges = linspace(minISI/scale,maxISI/scale,nbins);
    parfor i=1:n
        [isiMat(i,:), ~] = histcount(isi{i},edges); 
        isiMat(i,:)./max(isiMat(i,:))*100;
    end
end
function plotsub(isiMat,edges,scale)
    m = mean(isiMat);
    l = m-quantile(isiMat,0.25);
    u = quantile(isiMat,0.75) - m;
    edges = edges + (edges(2) - edges(1))/2;
    errorbar(edges*scale, m, l, u);
end
