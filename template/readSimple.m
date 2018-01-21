function [neuronProperty,tuningCurve] = readSimple(DIR,ntheta,n,ndperiod)
    filepath = [DIR, '/', 'frate.dat'];
    fid = fopen(filepath);
    tuningCurve.frate = fread(fid, [n, 2*ntheta],'double');
    tuningCurve.stifr = fread(fid, [n, 2*ntheta],'double');
    fclose(fid);
    
    filepath = [DIR, '/', 'cv.dat'];
    fid = fopen(filepath);
    Data = fread(fid, [n, 16],'double');
    neuronProperty.cv = Data(:,1);
    neuronProperty.sc = Data(:,2);
    neuronProperty.pkrate = Data(:,5);
    neuronProperty.prA = Data(:,7);
    neuronProperty.ei = Data(:,8);
    neuronProperty.br = Data(:,9);
    neuronProperty.cvNoBack = Data(:,10);
    % neuronProperty.oorate = Data(:,11);
    neuronProperty.priA = Data(:,12);
    neuronProperty.indpo = nearest(Data(:,13));
    neuronProperty.indoo = nearest(Data(:,14));
    neuronProperty.indpi = nearest(Data(:,15));
    neuronProperty.indoi = nearest(Data(:,16));
    fclose(fid);
    
    filepath = [DIR, '/', 'intra.dat'];
    fid = fopen(filepath);
    tuningCurve.gtot   = fread(fid, [n, 2*ntheta+1],'double');
                         fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.Veff   = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gLGN   = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gE     = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gI     = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gEn    = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gIn    = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.rate   = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.V      = fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gLGNsc = fread(fid, [n, 2*ntheta+1],'double');
                         fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gEsc   = fread(fid, [n, 2*ntheta+1],'double');
                         fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.Veffsc = fread(fid, [n, 2*ntheta+1],'double');
                         fread(fid, [n, 2*ntheta+1],'double');
    tuningCurve.gP     = fread(fid, [n, 2*ntheta+1],'double');
    fclose(fid);
end

