function [indpi,Veff,Veffstd,Veffsc,rate,allCycle] = readVsFR(DIR,ntheta,n,ndperiod)

filepath = [DIR, '/', 'cv.dat'];
fid = fopen(filepath);
fread(fid, [n, 11],'double');
priA = nearest(fread(fid, [n, 1],'double'));
dtheta = 180/ntheta;
priA = round(priA/pi*180/dtheta)+1;
fclose(fid);

filepath = [DIR, '/', 'intra.dat'];
fid = fopen(filepath);

fread(fid, [n,2*2*ntheta+1],'double');
Veff = fread(fid, [n, ntheta],'double');
fread(fid, [n,ntheta+1],'double');

fread(fid, [n,5*2*ntheta+1],'double');
rate = fread(fid,[n,ntheta],'double');
fread(fid,[n,ntheta+1],'double');

fread(fid, [n,5*2*ntheta+1],'double');
Veffsc = fread(fid, [n, ntheta],'double');

fclose(fid);

filepath = [DIR, '/', 'intra2.dat'];
fid = fopen(filepath);
fread(fid, [n,2*2*ntheta+1],'double');
Veff2 = fread(fid, [n, ntheta],'double');
fclose(fid);

Veffstd = sqrt(Veff2-Veff.^2);

filepath = [DIR, '/', 'cycles.dat'];
fid = fopen(filepath);
csp = zeros(n,ndperiod,ntheta);cgl=csp;cge=csp;cgi=csp;cgen=csp;cgin=csp;cvm=csp;cvs=csp;cal=csp;cbe=csp;cgp=csp;
cgl2=csp;cge2=csp;cgi2=csp;cgen2=csp;cgin2=csp;cvm2=csp;cvs2=csp;cal2=csp;cbe2=csp;cgp2=csp;
for i=1:ntheta
    csp(:,:,i) = fread(fid,[n, ndperiod],'double');
    cgl(:,:,i) = fread(fid,[n, ndperiod],'double');
    cge(:,:,i) = fread(fid,[n, ndperiod],'double');
    cgi(:,:,i) = fread(fid,[n, ndperiod],'double');
    cgen(:,:,i) = fread(fid,[n, ndperiod],'double');
    cgin(:,:,i) = fread(fid,[n, ndperiod],'double');
    cvm(:,:,i) = fread(fid,[n, ndperiod],'double');
    cal(:,:,i) = fread(fid,[n, ndperiod],'double');
    cbe(:,:,i) = fread(fid,[n, ndperiod],'double');
    cvs(:,:,i) = fread(fid,[n, ndperiod],'double');
    cgp(:,:,i) = fread(fid,[n, ndperiod],'double');
    fread(fid,[n, 6*ndperiod],'double');
end
for i=1:ntheta
	cgl2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cge2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cgi2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cgen2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cgin2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cvm2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cal2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cbe2(:,:,i) = fread(fid,[n, ndperiod],'double');
	cvs2(:,:,i) = fread(fid,[n, ndperiod],'double');
    cgp2(:,:,i) = fread(fid,[n, ndperiod],'double');
    fread(fid,[n, 6*ndperiod],'double');
end
fclose(fid);

allCycle.gtotstd = sqrt(cal2-cal.*cal);
allCycle.Itotstd = sqrt(cbe2-cbe.*cbe);
allCycle.VsSTD = sqrt(cvs2-cvs.*cvs);
allCycle.gLGNstd = sqrt(cgl2-cgl.*cgl);
allCycle.gEstd = sqrt(cge2-cge.*cge);
allCycle.gIstd = sqrt(cgi2-cgi.*cgi);
allCycle.gEnstd = sqrt(cgen2-cgen.*cgen);
allCycle.gInstd = sqrt(cgin2-cgin.*cgin);
allCycle.gPstd = sqrt(cgp2-cgp.*cgp);
allCycle.Vstd = sqrt(cvm2-cvm.*cvm);
allCycle.spikes = csp;
allCycle.gtot = cal;
allCycle.Itot = cbe;
allCycle.Vs = cvs;
allCycle.gLGN = cgl;
allCycle.gE = cge;
allCycle.gI = cgi;
allCycle.gEn = cgen;
allCycle.gIn = cgin;
allCycle.gP = cgp;
allCycle.V = cvm;
end
