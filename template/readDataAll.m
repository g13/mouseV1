function [neuronProperty,tuningCurve,allCycle] = readDataAll(DIR,ntheta,n,ndperiod)

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
tuningCurve.gtot = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.Itot = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.Veff = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gLGN = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gE = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gI = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gEn = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gIn = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.rate = fread(fid,[n,2*ntheta+1],'double');
tuningCurve.V = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gLGNsc = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gLGNsc2 = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gEsc = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gEsc2 = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.Veffsc = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.Veffsc2 = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.gP = fread(fid, [n, 2*ntheta+1],'double');
fclose(fid);

filepath = [DIR, '/', 'curr.dat'];
fid = fopen(filepath);
tuningCurve.cLGN = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.cLGN_F1 = zeros(n,ntheta);
tuningCurve.cE = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.cE_F1 = zeros(n,ntheta);
tuningCurve.cI = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.cEn = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.cIn = fread(fid, [n, 2*ntheta+1],'double');
cLGN2 = fread(fid, [n, 2*ntheta+1],'double');
cE2 = fread(fid, [n, 2*ntheta+1],'double');
cI2 = fread(fid, [n, 2*ntheta+1],'double');
cEn2 = fread(fid, [n, 2*ntheta+1],'double');
cIn2 = fread(fid, [n, 2*ntheta+1],'double');
tuningCurve.cP = fread(fid, [n, 2*ntheta+1],'double');
cP2 = fread(fid, [n, 2*ntheta+1],'double');
fclose(fid);
tuningCurve.cLGNstd = sqrt(cLGN2-tuningCurve.cLGN.^2);
tuningCurve.cEstd = sqrt(cE2-tuningCurve.cE.^2);
tuningCurve.cIstd = sqrt(cI2-tuningCurve.cI.^2);
tuningCurve.cEnstd = sqrt(cEn2-tuningCurve.cEn.^2);
tuningCurve.cInstd = sqrt(cIn2-tuningCurve.cIn.^2);
tuningCurve.cPstd = sqrt(cP2-tuningCurve.cP.^2);

filepath = [DIR, '/', 'intra2.dat'];
fid = fopen(filepath);
gtot2 = fread(fid, [n, 2*ntheta+1],'double');
Itot2 = fread(fid, [n, 2*ntheta+1],'double');
Veff2 = fread(fid, [n, 2*ntheta+1],'double');
gLGN2 = fread(fid, [n, 2*ntheta+1],'double');
gE2 = fread(fid, [n, 2*ntheta+1],'double');
gI2 = fread(fid, [n, 2*ntheta+1],'double');
gEn2 = fread(fid, [n, 2*ntheta+1],'double');
gIn2 = fread(fid, [n, 2*ntheta+1],'double');
rate2 = fread(fid,[n,2*ntheta+1],'double');
V2 = fread(fid, [n, 2*ntheta+1],'double');
gP2 = fread(fid, [n, 2*ntheta+1],'double');

fclose(fid);
tuningCurve.gtotstd = sqrt(gtot2-tuningCurve.gtot.^2);
tuningCurve.Itotstd = sqrt(Itot2-tuningCurve.Itot.^2);
tuningCurve.Veffstd = sqrt(Veff2-tuningCurve.Veff.^2);
tuningCurve.gLGNstd = sqrt(gLGN2-tuningCurve.gLGN.^2);
tuningCurve.gEstd = sqrt(gE2-tuningCurve.gE.^2);
tuningCurve.gIstd = sqrt(gI2-tuningCurve.gI.^2);
tuningCurve.gEnstd = sqrt(gEn2-tuningCurve.gEn.^2);
tuningCurve.gInstd = sqrt(gIn2-tuningCurve.gIn.^2);
tuningCurve.Vstd = sqrt(V2-tuningCurve.V.^2);
tuningCurve.gPstd = sqrt(gP2-tuningCurve.gP.^2);

filepath = [DIR, '/', 'cycles.dat'];
fid = fopen(filepath);
csp = zeros(n,ndperiod,ntheta);cgl=csp;cge=csp;cgi=csp;cgen=csp;cgin=csp;cvm=csp;cvs=csp;cal=csp;cbe=csp;cgp=csp;
cgl2=csp;cge2=csp;cgi2=csp;cgen2=csp;cgin2=csp;cvm2=csp;cvs2=csp;cal2=csp;cbe2=csp;cgp2=csp;
clgn=csp;cexc=csp;cinh=csp;cfexc=csp;cfinh=csp;cpota=csp;
clgn2=csp;cexc2=csp;cinh2=csp;cfexc2=csp;cfinh2=csp;cpota2=csp;
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
    clgn(:,:,i) = fread(fid,[n, ndperiod],'double');
    cexc(:,:,i) = fread(fid,[n, ndperiod],'double');
    cinh(:,:,i) = fread(fid,[n, ndperiod],'double');
    cfexc(:,:,i) = fread(fid,[n, ndperiod],'double');
    cfinh(:,:,i) = fread(fid,[n, ndperiod],'double');
    cpota(:,:,i) = fread(fid,[n, ndperiod],'double');
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
    clgn2(:,:,i) = fread(fid,[n, ndperiod],'double');
    cexc2(:,:,i) = fread(fid,[n, ndperiod],'double');
    cinh2(:,:,i) = fread(fid,[n, ndperiod],'double');
    cfexc2(:,:,i) = fread(fid,[n, ndperiod],'double');
    cfinh2(:,:,i) = fread(fid,[n, ndperiod],'double');
    cpota2(:,:,i) = fread(fid,[n, ndperiod],'double');
end
fclose(fid);

allCycle.Vstd = sqrt(cvm2-cvm.*cvm);
allCycle.gtotstd = sqrt(cal2-cal.*cal);
allCycle.Itotstd = sqrt(cbe2-cbe.*cbe);
allCycle.Veffstd = sqrt(cvs2-cvs.*cvs);
allCycle.gLGNstd = sqrt(cgl2-cgl.*cgl);
allCycle.gEstd = sqrt(cge2-cge.*cge);
allCycle.gIstd = sqrt(cgi2-cgi.*cgi);
allCycle.gEnstd = sqrt(cgen2-cgen.*cgen);
allCycle.gInstd = sqrt(cgin2-cgin.*cgin);
allCycle.gPstd = sqrt(cgp2-cgp.*cgp);
allCycle.Vstd = sqrt(cvm2-cvm.*cvm);
allCycle.cLGNstd= sqrt(clgn2-clgn.*clgn);
allCycle.cEstd= sqrt(cexc2-cexc.*cexc);
allCycle.cIstd= sqrt(cinh2-cinh.*cinh);
allCycle.cEnstd=sqrt(cfexc2-cfexc.*cfexc);
allCycle.cInstd=sqrt(cfinh2-cfinh.*cfinh);
allCycle.cPstd= sqrt(cpota2-cpota.*cpota);
allCycle.spikes = csp;
allCycle.gtot = cal;
allCycle.Itot = cbe;
allCycle.Veff = cvs;
allCycle.gLGN = cgl;
allCycle.gE = cge;
allCycle.gI = cgi;
allCycle.gEn = cgen;
allCycle.gIn = cgin;
allCycle.gP = cgp;
allCycle.V = cvm;
allCycle.cLGN=clgn;
allCycle.cE=cexc;
allCycle.cI=cinh;
allCycle.cEn=cfexc;
allCycle.cIn=cfinh;
allCycle.cP=cpota;
end
