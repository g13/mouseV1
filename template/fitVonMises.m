function [r0, rmax, smax, D, adjrs] = fitVonMises(frate,prA,levels,n,theme,ntheta,threads,DIR)
    rmax = zeros(levels,n);
    r0 = rmax;
    smax = rmax;
    D = rmax;
    adjrs = rmax;
    dtheta = pi/ntheta;
    hWF=sqrt(log(2));  %halfWidthFactor
    pool = gcp('nocreate');
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

    Gauss1 = fittype(...
		     @(smax,D,rmax,r0,x) r0+rmax*exp((cos(2*(x-smax))-1)./D),...
			 'independent',{'x'},'coefficients',{'smax','D','rmax','r0'});
    for ii = 1:levels
		disp(['fitting ',num2str(ii)]);
        parfor (j = 1:n,threads)
            ipA = round(prA(j,ii)/dtheta)+1;
            if ipA < ntheta/2+1
                shift = [(ipA+ntheta/2):ntheta,1:(ipA+ntheta/2)];
            else 
                shift = [(ipA-ntheta/2):ntheta,1:(ipA-ntheta/2)];
            end
            thetas = ((ipA-ntheta/2-1):(ipA+ntheta/2-1))*dtheta;
            smaxRight = thetas(ntheta*3/4);
            smaxLeft = thetas(ntheta/4);
            fr = frate(j,shift,ii);
            rmax0 = max(fr)-min(fr);
            stdfr = std(fr);
            [fitted, gof] = fit(thetas', fr',Gauss1,...
                'Lower',[smaxLeft,				    0,	    0,        0],...
				'Upper',[smaxRight,					inf,    1.5*max(fr), max(fr)],...
                'Start',[(thetas(1)+thetas(end))/2, 1e4,	mean(fr), mean(fr)]);

            r0(ii,j) = fitted.r0;
            rmax(ii,j) = fitted.rmax;
            smax(ii,j) = fitted.smax;
            D(ii,j) = fitted.D;
            adjrs(ii,j) = gof.adjrsquare;
        end
    end
    save([DIR,'/',theme,'-fitted.mat'], 'r0','rmax','smax','D','adjrs');
end
