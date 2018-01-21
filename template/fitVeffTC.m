function [vmax, smax, sigfsq2, lifted] = fitVeffTC(Veff,prA,levels,n,ntheta,threads)
    vmax = zeros(1,n)-1;
    smax = vmax;
    sigfsq2 = vmax;
    dtheta = 180/ntheta;
    for ii = 1:levels
		disp(['fitting Veff for contrast ',num2str(ii)]);

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
            	vmax(j) = 0;
            	smax(j) = (thetas(1)+thetas(end))/2;
            	sigfsq2(j) = 180;
            	lifted(j) = 0;
			else
            	if (max(vs)-mean(vs))<2*std(vs)
            	    vmax0 = 0;
            	else
            	    vmax0 = max(vs)-min(vs);
            	end
				lifted0 = min(0,min(vs));
            	liftedGauss1 = fittype(...
					@(smax,sigfsq2,vmax,lifted,x) lifted + vmax*exp(-((x-smax)./sigfsq2).^2),...
					'independent',{'x'},'coefficients',{'smax','sigfsq2','vmax','lifted'});
            	fitted = fit(thetas', vs',liftedGauss1,...
            	    'Lower',[smaxLeft,					1e-8,	max(vs)-mean(vs),								lifted0],...
					'Upper',[smaxRight,					Inf,	max(max(vs)-mean(vs)+2e-7,2*(max(vs)-min(vs))),	mean(vs)],...
            	    'Start',[(thetas(1)+thetas(end))/2,	40,		vmax0,											lifted0]);
            	vmax(j) = fitted.vmax;
            	smax(j) = fitted.smax;
            	sigfsq2(j) = fitted.sigfsq2;
            	lifted(j) = fitted.lifted;
			end
        end
    end
end
