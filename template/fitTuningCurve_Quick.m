function [rmax, smax, sigfsq2, lifted] = fitTuningCurve_Quick(frate,br,prA,contrast,n,theme,ntheta,threads,substract,DIR)
	levels = length(contrast);
    rmax = zeros(1,n)-1;
    smax = rmax;
    sigfsq2 = rmax;
    dtheta = 180/ntheta;
    for ii = 1:levels
		if substract
			disp(['fitting subtracted ',num2str(ii)]);
		else
			disp(['fitting ',num2str(ii)]);
		end
		i = contrast(ii);

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
			if (sum(fr) == 0)
            	rmax(j) = 0;
            	smax(j) = (thetas(1)+thetas(end))/2;
            	sigfsq2(j) = 180;
            	lifted(j) = 0;
			else
				if substract
            	    fr = fr-br(j,ii);
            	end
            	if (max(fr)-mean(fr))<2*std(fr)
            	    rmax0 = 0;
            	else
            	    rmax0 = max(fr)-min(fr);
            	end
				lifted0 = min(0,min(fr));
            	liftedGauss1 = fittype(...
					@(smax,sigfsq2,rmax,lifted,x) lifted + rmax*exp(-((x-smax)./sigfsq2).^2),...
					'independent',{'x'},'coefficients',{'smax','sigfsq2','rmax','lifted'});
            	fitted = fit(thetas', fr',liftedGauss1,...
            	    'Lower',[smaxLeft,					1e-8,	max(fr)-mean(fr),								lifted0],...
					'Upper',[smaxRight,					Inf,	max(max(fr)-mean(fr)+2e-7,2*(max(fr)-min(fr))),	mean(fr)],...
            	    'Start',[(thetas(1)+thetas(end))/2,	40,		rmax0,											lifted0]);
            	rmax(j) = fitted.rmax;
            	smax(j) = fitted.smax;
            	sigfsq2(j) = fitted.sigfsq2;
            	lifted(j) = fitted.lifted;
			end
        end

		if substract
    	    save([DIR,'/',theme,'-',num2str(i),'fitted_substracted.mat'], 'rmax','smax','sigfsq2','lifted');
    	else
    	    save([DIR,'/',theme,'-',num2str(i),'fitted.mat'], 'rmax','smax','sigfsq2','lifted');
    	end
    end
end
