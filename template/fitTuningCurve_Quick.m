function [rmax, smax, sigfsq2] = fitTuningCurve_Quick(frate,prA,levels,n,theme,ntheta,threads,DIR)
%function [rmax, smax, sigfsq2, lifted] = fitTuningCurve_Quick(frate,cv,prA,levels,n,theme,ntheta,threads,DIR)
    rmax = zeros(levels,n);
    smax = rmax;
    sigfsq2 = rmax;
    dtheta = 180/ntheta;
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

    for ii = 1:levels
		disp(['fitting ',num2str(ii)]);
        parfor (j = 1:n,threads)
            ipA = round(prA(j,ii)*180/pi/dtheta)+1;
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
            Gauss1 = fittype(...
				@(smax,sigfsq2,rmax,x) rmax*exp(-((x-smax)./sigfsq2).^2),...
				'independent',{'x'},'coefficients',{'smax','sigfsq2','rmax'});
            stdfr = std(fr);
            fitted = fit(thetas', fr',Gauss1,...
                'Lower',[smaxLeft,				    dtheta/(2*hWF),	0],...
				'Upper',[smaxRight,					inf,	        1.5*max(fr)],...
                'Start',[(thetas(1)+thetas(end))/2,	stdfr*sqrt(2),	max(fr)]);

            rmax(ii,j) = fitted.rmax;
            smax(ii,j) = fitted.smax;
            sigfsq2(ii,j) = fitted.sigfsq2;
        end
    end
    save([DIR,'/',theme,'-fitted.mat'], 'rmax','smax','sigfsq2');
end
