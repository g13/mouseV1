function [vmax, smax, sigfsq2, lifted] = fitVeffTC(Veff,prA,levels,n,ntheta,threads)
    vmax = zeros(n,levels);
    smax = vmax;
    sigfsq2 = vmax;
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
		disp(['fitting Veff for contrast ',num2str(ii)]);
		disp(['expecting ', num2str(sum(sum(Veff(:,:,ii),2)==0)),' silences']);

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
            vs = Veff(j,shift,ii);
			if (sum(vs) == 0)
            	vmax(j,ii) = 0;
            	smax(j,ii) = (thetas(1)+thetas(end))/2;
            	sigfsq2(j,ii) = inf;
            	lifted(j,ii) = 0;
			else
            	vmax0 = max(vs)-min(vs);
            	liftedGauss1 = fittype(...
					@(smax,sigfsq2,vmax,lifted,x) lifted + vmax*exp(-((x-smax)./sigfsq2).^2),...
					'independent',{'x'},'coefficients',{'smax','sigfsq2','vmax','lifted'});
                if max(vs)-mean(vs) > 1.45*std(vs)
            	    fitted = fit(thetas', vs',liftedGauss1,...
            	        'Lower',[smaxLeft,					dtheta/hWF,	max(vs)-mean(vs),		0],...
				    	'Upper',[smaxRight,					inf,	    2*(max(vs)-min(vs)),	mean(vs)],...
            	        'Start',[(thetas(1)+thetas(end))/2,	45/hWF,		vmax0,					min(vs)]);
                else
            	    fitted = fit(thetas', vs',liftedGauss1,...
            	        'Lower',[smaxLeft,					45/hWF,	  max(vs)-mean(vs),		   0],...
				    	'Upper',[smaxRight,					inf,	    1.5*(max(vs)-min(vs)),	mean(vs)],...
            	        'Start',[(thetas(1)+thetas(end))/2,	1e10,		0,					    mean(vs)]);

                end
            	vmax(j,ii) = fitted.vmax;
            	smax(j,ii) = fitted.smax;
            	sigfsq2(j,ii) = fitted.sigfsq2;
            	lifted(j,ii) = fitted.lifted;
			end
        end
    end
end
