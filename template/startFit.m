function startFit(nv1,theme,contrast,ntheta,threads,cleanclose)
    nv1=10800;
    ndperiod = 25;
	disp([num2str(length(contrast)),' contrasts']);
	disp([num2str(ntheta),' drifting orientations']);
	disp(['using ', num2str(threads), ' threads']);
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
    nc = length(contrast);
    br = zeros(nv1,nc);
    prA = zeros(nv1,nc);
    frate = zeros(nv1,2*ntheta,nc);
	for ii=1:nc
		i = contrast(ii);
	    eps = 125*2^(i-1);
		DIR = [theme,'/',num2str(eps,'%0.4d')];
		disp(DIR);
        filepath = [DIR, '/', 'frate.dat'];
        fid = fopen(filepath);
        frate(:,:,ii) = fread(fid, [nv1, 2*ntheta],'double');
        fclose(fid);

        filepath = [DIR, '/', 'cv.dat'];
        fid = fopen(filepath);
        Data = fread(fid, [nv1, 16],'double');
        prA(:,ii)=Data(:,7);
        br(:,ii)=Data(:,9);
        fclose(fid);
		%[nP(ii),tC(ii),~] = readDataAll(DIR,ntheta,nv1,ndperiod);
	end
	fitTuningCurve_Quick(frate,br,prA,contrast,nv1,theme,ntheta,threads,false,theme);
	%fitTuningCurve(tC,nP,contrast,nv1,theme,ntheta,threads,false,theme);
	%fitTuningCurve(tC,nP,contrast,nv1,theme,ntheta,threads,true,[theme,'/processed']);
	if cleanclose && threads > 1
		delete(pool);
	end
end
