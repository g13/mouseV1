%% generate connection matrix
function generateConMat(suffix,theme,preMatProfile,reciprocal,eSpecific,seed,plotout,format,scType,heterogeneous,outputName,eiSpecific,ieSpecific,coMatfile,logP)
%!!!!!!!!!!!!!!!!!! for scType == 'n' decide connection by relative
%cortical strength vs lgn strength and portion in total excitation
if ~exist('conMatFigures','dir')
    mkdir('conMatFigures');
end
assert(exist('conMatFigures'));
pPosition = [18, 180, 1200, 900];
if nargin < 15
    disp('no logNormal, return');
    return
    if nargin <14
        coMatfile ='';
        if nargin < 13
            ieSpecific = false;
            if nargin < 12
                eiSpecific = false;
                if nargin <11
                    outputName = '';
                    if nargin < 10
                        heterogeneous = '';
                        if nargin < 9
                            scType = 'n';
                            if nargin < 8
                                format = '';
                                if nargin < 7
                                    plotout = false;
                                    if nargin < 6
                                        seed = 0;
                                        if nargin < 5
                                            eSpecific = 'none';
                                            if nargin < 4
                                                reciprocal = 0.0;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
	if ~isempty(format)
		if strcmp(format,'psc2')
			printDriver = ['-de',format];
			format = 'eps';
		else
			printDriver = ['-d',format];
		end
		dpi = '-r150';
	end
FontSize = 20;
set(0,'DefaultAxesFontSize',FontSize)
switch preMatProfile
    case 'uniform'
        mcase = 0;
        bounde = sqrt(2);
        boundi = sqrt(2);
    case 'gauss'
        mcase = 1;
        bounde = 0;
        boundi = sqrt(2);
    otherwise
        disp(['no function is defined for ', preMatProfile]);
        return
end
switch heterogeneous
    case 'logNormal'
        q.profile = @logNormalProfile;
    case 'None'
        q.profile = function_handle.empty;
    otherwise
        disp(['no function is defined for ', heterogeneous]);
end
matName = [suffix,'-',theme,'-',preMatProfile,'-',eSpecific,...
            '-r',num2str(reciprocal*100,'%d'),'-s',num2str(seed),scType,'-',heterogeneous];
load([suffix,'.mat']);
% if ~strcmp(eSpecific, '')
% 	load([suffix,'.mat']);
%     etheta = reshape(etheta,[p.nv1e,1]);
% end
if strcmp(eSpecific,'coGauss')
    eeSigCoeff = 0.8;
    if eiSpecific
        eiSigCoeff = 1.00;
    end
    if ieSpecific
        ieSigCoeff = 0.50;
    end
    %coMatfile = ['coMa-60x60-1xu-nmLDoO-ff3-s911.mat'];
    load(coMatfile);
    size(coMat)
end
if logP.spread
    logP.normDistance = p.enormDistance;
end
clear p.eSubgion p.esigma p.ePeakDistance p.isigma p.iPeakDistance p.iAspectRatio p.eAspectRatio p.iSubregion;
elx = (p.ev1x) * p.dE;
ely = (p.ev1y) * p.dE;
ilx = (p.iv1x) * p.dI;
ily = (p.iv1y) * p.dI;

% r = Half-Width-Half-Maximum

e.raxn = 100;
e.rden = 75;
i.raxn = 80;
i.rden = 50;
e.sige = 1;
e.sigi = 1;
i.sige = 1;
i.sigi = 1;
%i.raxn = 90;
%i.rden = 50;

% e.raxn = 80;
% e.rden = 80;
% i.raxn = 80;
% i.rden = 80;

e.r = p.nv1e/p.nv1;
i.r = 1-e.r;

e.probe = 0.1765;
e.probi = 0.1765;
i.probe = 0.6;
i.probi = 0.9;
% i.probe = 0.1765;
% i.probi = 0.706;

%i.probe = 0.40;
%i.probi = 1.0;

if reciprocal > 0.0
    rcase = mcase;
end

self = false; % self connection
q.exact = true; % homogeneouse incoming connections.

[nee0, nei0, nie0, nii0] = getNeff(e,i,'2_3');
nex = 5;
disp(['nee0 = ',num2str(nee0)]);
% statistics correction of rand num generation deviation
%switch preMatProfile
%    case 'uniform'
%		nee = round(nee0 + 7);
%        nei = round(nei0 - 9);
%        nie = round(nie0 - 15);
%        nii = round(nii0 + 47);
%    case 'gauss'
		nee = round(nee0 + 10);
        nei = round(nei0 + 2);
        nie = round(nie0 + 32);
        nii = round(nii0 + 26);
%    otherwise
%end

nLGNmean = mean(nLGN(1:p.nv1e,1));

% strength ratio LGN over CC
sRatio = (11.0/nLGNmean) / (0.45);
% rate estimation ratio LGN over CC
rRatio = 50/6;

ratio = rRatio * sRatio;
% portion of excitation
%portion = 3/7;
portion = 407/800;

nEEmean = nLGNmean*ratio/portion;
totalCon = nEEmean/ratio + nLGNmean;
disp('prescibed pre:');

switch scType 
	case 'n'
        correction = 0;
        %compressedLGN = nLGN(1:p.nv1e,1);
        compressedLGN = ones(p.nv1e,1)*nLGNmean;
        compressedLGN = nLGNmean + 0.5*(compressedLGN-nLGNmean);
         nee = round(ratio*(totalCon  - compressedLGN)) + correction;
       %nee = round(ratio*(totalCon - mean(nLGN(1:p.nv1e,1)))) + correction;
        nee(nee<0) = 0;
        q.rep = 4;
        q.repRetro = q.rep;
        nei = round(nEEmean * e.probi/(4*e.probe));
        nie = round(nie0);
        nii = round(nii0);
        disp(['maxEE=',num2str(max(nee)),', minEE=',num2str(min(nee))]);
        disp(['maxLGN=',num2str(max(nLGN(1:p.nv1e,1))),', minLGN=',num2str(min(nLGN(1:p.nv1e,1)))]);
        
        disp(['EE:',num2str(round(nEEmean)),', EI:',num2str(nei)]);
        disp(['IE:',num2str(nie),', II:',num2str(nii)]);
	case 's'
        q.rep = 1;
        q.repRetro = q.rep;
        disp(['EE:',num2str(round(nee0)),', EI:',num2str(round(nei0))]);
        disp(['IE:',num2str(round(nie0)),', II:',num2str(round(nii0))]);
	otherwise
end

pex = ones(p.ev1y,1)*((1:p.ev1x)-0.5).*p.dE;
pex = reshape(pex,[p.nv1e,1]); % id goes along y first (column major)
pey = ((1:p.ev1y)-0.5)'*ones(1,p.ev1x).*p.dE;
pey = reshape(pey,[p.nv1e,1]);

pix = ones(p.iv1y,1) * (((1:p.iv1x)-0.5).*p.dI);
pix = reshape(pix,[p.nv1i,1]);
piy = (((1:p.iv1y)-0.5)'.*p.dI) * ones(1,p.iv1x);
piy = reshape(piy,[p.nv1i,1]);

% We choose to find presynaptic neurons
% fix one neuron then try to find e or i, no matter pre or post 
ae = p.dE*p.dE;
ai = p.dI*p.dI;
tic;
%% initial connection matrix
% (i,j) => i innervate j
m = zeros(p.nv1,p.nv1,'int8');

% EE
if seed ~= 0, seed = seed+1; end
% pre E
src.x = pex;
src.y = pey;
src.n = p.nv1e;
src.a = ae;
src.neff = nee;
src.lx = elx;
src.ly = ely;
src.raxn = e.raxn;
% post E
tar.x = pex;
tar.y = pey;
tar.n = p.nv1e;
tar.rden = e.rden;
tar.bound = bounde;
tar.sig = e.sige;
q.self = self;
switch eSpecific
%    case 'thetaGauss'
%        tar.sigtheta = 40/(2*sqrt(2*log(2))); % in degree
%		src.theta = etheta;
%		tar.theta = etheta;
%        q.self = false;
%		switch preMatProfile
%		    case 'uniform'
%                
%		    case 'gauss'
%				tar.bound = 0;
%		    otherwise
%		        disp(['thetaGauss function is not defined for ', preMatProfile]);
%		        return
%		end
%    case 'thetaUniform'
%        src.theta = etheta;
%		tar.theta = etheta;
%        switch preMatProfile
%		    case 'uniform'
%                
%		    case 'gauss'
%				tar.bound = 0;
%		    otherwise
%		        disp(['thetaGauss function is not defined for ', preMatProfile]);
%		        return
%        end
    case 'coGauss'
        tar.sigCoeff = eeSigCoeff;
        q.specificMat = 1-coMat(1:p.nv1e,1:p.nv1e);
        switch preMatProfile
		    case 'uniform'
                q.case = 2;    
            case 'gauss'
                q.case = 3;
            otherwise
		        disp(['coGauss function is not defined for ', preMatProfile]);
        end
    case 'coUniform'
        switch preMatProfile
		    case 'uniform'
                q.case = 4;    
            case 'gauss'
                q.case = 5;
            otherwise
		        disp(['coUniform function is not defined for ', preMatProfile]);
        end
    case 'none'
        q.case = mcase;
end
tic;
if logP.spread
    [mee,h] = subMat(src,tar,q,logP); 
    if ~isempty(format)
        set(h, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(h,['EPSPoverND','-',matName,'.',format],printDriver,dpi);
    end
else
    mee = subMat(src,tar,q,logP); 
end
disp('EE finished');
toc;
m(1:p.nv1e,1:p.nv1e) = mee;

if isempty(outputName)
    save(['mee-',matName,'.mat'],'mee');
else
    save(['mee-',outputName,'.mat'],'mee');
end
q.profile = function_handle.empty();
%EI

% pre I
src.x = pix;
src.y = piy;
src.n = p.nv1i;
src.a = ai;
src.neff = nei;
src.lx = ilx;
src.ly = ily;
src.raxn = i.raxn;
% post E
tar.x = pex;
tar.y = pey;
tar.n = p.nv1e;
tar.rden = e.rden;
tar.bound = bounde;
tar.sig = e.sigi;
q.case = mcase;
q.self = false;

casetmp = q.case;
% tar.bound = sqrt(2);
% q.case = 2;
if eiSpecific
    switch eSpecific
        case 'coGauss'
            tar.sigCoeff = eiSigCoeff;
            q.specificMat = 1-coMat(p.nv1e+(1:p.nv1i),1:p.nv1e);
            switch preMatProfile
                case 'uniform'
                    q.case = 2;    
                case 'gauss'
                    q.case = 3;
                otherwise
                    disp(['coGauss function is not defined for ', preMatProfile]);
            end
        case 'coUniform'
            switch preMatProfile
	    	    case 'uniform'
                    q.case = 4;    
                case 'gauss'
                    q.case = 5;
                otherwise
	    	        disp(['coUniform function is not defined for ', preMatProfile]);
            end
        case 'none'
            q.case = mcase;
    end
end
mei = subMat(src,tar,q);
if eiSpecific
    disp('eiSpecific finished');
end
q.case = casetmp;

m(p.nv1e+1:p.nv1,1:p.nv1e) = mei;
%IE

% pre E
src.x = pex;
src.y = pey;
src.n = p.nv1e;
src.a = ae;
src.neff = nie;
src.lx = elx;
src.ly = ely;
src.raxn = e.raxn;
% post I
tar.x = pix;
tar.y = piy;
tar.n = p.nv1i;
tar.rden = i.rden;
tar.bound = boundi;
tar.sig = i.sige;
q.case = mcase;
q.self = false;

if reciprocal > 0.0
    q.rcase = rcase;
    [mie, ~] = retroMat_beta(src,tar,q,reciprocal,mei);
else
    if ieSpecific
        switch eSpecific
            case 'coGauss'
                tar.sigCoeff = ieSigCoeff;
                q.specificMat = 1-coMat(1:p.nv1e,p.nv1e+(1:p.nv1i));
                switch preMatProfile
                    case 'uniform'
                        q.case = 2;    
                    case 'gauss'
                        q.case = 3;
                    otherwise
                        disp(['coGauss function is not defined for ', preMatProfile]);
                end
            case 'coUniform'
                switch preMatProfile
	        	    case 'uniform'
                        q.case = 4;    
                    case 'gauss'
                        q.case = 5;
                    otherwise
	        	        disp(['coUniform function is not defined for ', preMatProfile]);
                end
            case 'none'
                q.case = mcase;
        end
    end
    mie = subMat(src,tar,q);
end
m(1:p.nv1e,p.nv1e+1:p.nv1) = mie;

if isempty(outputName)
    save(['mei-',matName,'.mat'],'mei');
else
    save(['mei-',outputName,'.mat'],'mei');
end

%II

% pre I
src.x = pix;
src.y = piy;
src.n = p.nv1i;
src.a = ai;
src.neff = nii;
src.lx = ilx;
src.ly = ily;
src.raxn = i.raxn;
% post I
tar.x = pix;
tar.y = piy;
tar.n = p.nv1i;
tar.rden = i.rden;
tar.bound = boundi;
tar.sig = i.sigi;
q.case = mcase;
q.self = self;

mii = subMat(src,tar,q);
m(p.nv1e+1:p.nv1,p.nv1e+1:p.nv1) = mii;
toc;
if isempty(outputName)
    fid = fopen(['conMat-',matName,'.mat'],'w');
else
    fid = fopen([outputName,'.mat'],'w');
end
fwrite(fid,m,'int8');
fclose(fid);
if strcmp(eSpecific,'coGauss')
    clear coMat;
end
if plotout

eg = [randi(p.nv1e,[nex-1,1]), randi(p.nv1i,[nex-1,1])];
eg(nex,1) = p.nv1e/2+p.ev1y/2;
eg(nex,2) = p.nv1i/2+p.iv1y/2;
%% plot estimated RF situation
load('logNormalProfile.mat');
frtlgn0 = 3.24;
frtcE0 = 3;
frtcE = [frtcE0,profiles(2:end)./max(profiles(2:end))*2*frtcE0];
frtcI = 12;
SI = 0.08*4;
SE = 0.15;
g0 = 10;


rfx = @(x0,a,b,theta,t,ra,rb) x0 + (ra.*a.*cos(t).*cos(theta)-rb.*b.*sin(t).*sin(theta));
rfy = @(y0,a,b,theta,t,ra,rb) y0 + (ra.*a.*cos(t).*sin(theta)+rb.*b.*sin(t).*cos(theta));
FWHM = 2*sqrt(2*log(2));

ttt = 0:0.1:2*pi;
for j=1:nex
    i = p.nv1e + eg(j,2);
    hExRF=figure;
    hold on
    subregion = length(nSubLGN{i});
    subpick = 1:subregion;
    if i <= p.nv1e
        theta = etheta(i);
		sigma = p.esigma(i,subpick);
    	sigmb = sigma.*p.eAspectRatio(i,subpick);
    else
        theta = itheta(i-p.nv1e);
		sigma = p.isigma(i-p.nv1e,subpick);
    	sigmb = sigma.*p.iAspectRatio(i-p.nv1e,subpick);
    end
    if theta >=pi/2
        gtheta = theta - pi/2;
    else
        gtheta = theta + pi/2;
    end
    peak = zeros(subregion,2);
    for ij = 1:subregion
        if nSubLGN{i}(ij) > 0
            peak(ij,:) = mean(LGNpos(v1Map{i,ij},:),1)*180/pi;
        else
            peak(ij,:) = mean(LGNpos(-v1Map{i,ij},:),1)*180/pi;
        end
    end
    
    ymax = -inf;
    ymin = inf;
    xmax = ymax;
    xmin = ymin;
    for ij = 1:subregion
        xxx2 = rfx(peak(ij,1)/180*pi,sigmb(ij),sigma(ij),gtheta,ttt,FWHM,FWHM)*180/pi;         
        yyy2 = rfy(peak(ij,2)/180*pi,sigmb(ij),sigma(ij),gtheta,ttt,FWHM,FWHM)*180/pi;
        ymin = min(min(yyy2),ymin);
        ymax = max(max(yyy2),ymax);
        xmin = min(min(xxx2),xmin);
        xmax = max(max(xxx2),xmax);
    end
    dxy = p.rlgn*180/pi;
    edge = 5;
    x = (xmin-edge*dxy):(dxy/5):(xmax+edge*dxy);
    y = (ymin-edge*dxy):(dxy/5):(ymax+edge*dxy);
    scale = 700;
    zmax = 0;
    Z_ff_on = zeros(length(y),length(x));
    Z_ff_off = Z_ff_on;
    for ij = 1:subregion
  
        s = g0/(nLGNmean*frtlgn0);%*p.si;
        frtlgn = 50;
       if nSubLGN{i}(ij) > 0
            [X, Y, Z] = multiLGNspatialKernel(LGNpos(v1Map{i,ij},1)*180/pi, LGNpos(v1Map{i,ij},2)*180/pi,x,y,nSubLGN{i}(ij),lgnStrength{i,ij},s/scale*frtlgn);
            Z_ff_on = Z_ff_on + Z;
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','m','LineWidth',2);
            newZmax = max(max(Z));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle','--','LineColor','m','LineWidth',2);
       else
            [X, Y, Z] = multiLGNspatialKernel(LGNpos(-v1Map{i,ij},1)*180/pi, LGNpos(-v1Map{i,ij},2)*180/pi,x,y,-nSubLGN{i}(ij),lgnStrength{i,ij},s/scale*frtlgn);
            Z_ff_off = Z_ff_off + Z;
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','g','LineWidth',2);
            newZmax = max(max(Z));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle','--','LineColor','g','LineWidth',2);
       end
    end
%     mesh(X,Y,Z_ff_on,'FaceColor','none','EdgeColor','r');
%     mesh(X,Y,Z_ff_off,'FaceColor','none','EdgeColor','b');
    zmax = max(max(max(Z_ff_on)),max(max(Z_ff_off)));
    Z_Eon = zeros(size(X));
    Z_Eoff = Z_Eon;
    Z_Ion = Z_Eon;
    Z_Ioff = Z_Ion;
%       scale = 1200;
%     pick = 1:p.nv1e;
    pick = 1:p.nv1;
    neighbors = find(m(pick,i)>0);
    for ik = 1:length(neighbors)
        k = neighbors(ik);
%         assert(m(k,i)==1);
        if k <= p.nv1e      
            if m(k,i)== 1
                s = SE/nLGN(k,1);
            else
                s = profiles(m(k,i))/nLGN(k,1)/p2s;
            end
            frtc = frtcE(m(k,i));
        else
            frtc = frtcI;
            s = SI/nLGN(k,1);
        end
        subregion = length(nSubLGN{k});
        for ij = 1:subregion
           if nSubLGN{k}(ij) > 0
                [X, Y, Z] = multiLGNspatialKernel(LGNpos(v1Map{k,ij},1)*180/pi, LGNpos(v1Map{k,ij},2)*180/pi,x,y,nSubLGN{k}(ij),lgnStrength{k,ij},s/scale*frtc);
                if k <= p.nv1e
                    Z_Eon = Z_Eon + Z;
                else
                    Z_Ion = Z_Ion - Z;
                end
                %contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','m','LineWidth',2);
           else
                [X, Y, Z] = multiLGNspatialKernel(LGNpos(-v1Map{k,ij},1)*180/pi, LGNpos(-v1Map{k,ij},2)*180/pi,x,y,-nSubLGN{k}(ij),lgnStrength{k,ij},s/scale*frtc);
                if k <= p.nv1e
                    Z_Eoff = Z_Eoff + Z;
                else
                    Z_Ioff = Z_Ioff - Z;
                end
                %contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','g','LineWidth',2);
           end
%            mesh(X,Y,Z,'FaceColor','none');
        end
    end
    contour(X,Y,Z_Eon + Z_ff_on,[0,0],'LineStyle','--','LineColor','r','LineWidth',2);
    contour(X,Y,Z_Eoff + Z_ff_off,[0,0],'LineStyle','--','LineColor','b','LineWidth',2);
    Z_on = Z_Eon + Z_Ion + Z_ff_on;
    Z_off = Z_Eoff + Z_Ioff + Z_ff_off;
    zcmax = max(max(max(Z_on)),max(max(Z_off)));
    
    mesh(X,Y,Z_on-Z_off,'FaceColor','none','EdgeColor','k');
%     mesh(X,Y,Z_Eon+Z_ff_on,'FaceColor','none','EdgeColor','r');
%     mesh(X,Y,Z_Eoff+Z_ff_off,'FaceColor','none','EdgeColor','b');
%     mesh(X,Y,Z_Ion,'FaceColor','none','EdgeColor','r');
%     mesh(X,Y,Z_Ioff,'FaceColor','none','EdgeColor','b');
    
    [maxtmp, indi] = max(abs(Z_on-Z_off));
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'*k');
    title(['LGN RF max = ',num2str(zmax,'%3.1f'),', cortical RF max = ', num2str(zcmax,'%3.1f'),', max\delta=',num2str(deltamax,'%3.1f')]);
    
    [maxtmp, indi] = max(Z_on);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'*r');
    
    [maxtmp, indi] = max(Z_off);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'*b');
    
    [maxtmp, indi] = max(Z_Eon);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'^r');
    
    [maxtmp, indi] = max(Z_Eoff);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'^b');
    
    [maxtmp, indi] = max(Z_Ion);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'vr');
    
    [maxtmp, indi] = max(Z_Ioff);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'vb');
    
    [maxtmp, indi] = max(Z_ff_on);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'*m');
    
    [maxtmp, indi] = max(Z_ff_off);
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'*g');
    
    ylim([y(1),y(end)]);
    xlim([x(1),x(end)]);
    axis equal
    xlabel('degree');
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hExRF,['conMatFigures/egRF',num2str(j),'-',matName,'.',format],printDriver,dpi);
        saveas(hExRF,['conMatFigures/egRF',num2str(j),'-',matName,'.fig']);
    end
end
%%
for j = 1:nex
    hEx = figure;
    i = eg(j,1);
    subplot(2,2,1);
    hold on
    % E
%     x0 = p.ev1x/2;
%     y0 = p.ev1y/2;
%     i = (x0-1)*p.ev1y+y0;
    y0 = pey(i);
    x0 = pex(i);   
    % E<-E
    pick = m(1:p.nv1e,i)>0;
    x = pex(pick);
    y = pey(pick);
    npick = sum(pick);
    dx = zeros(npick,1);
    dy = zeros(npick,1);
    % enforce periodic boundary condition
            %x
            sel = abs(x-x0) > elx/2;
            dx(sel) = 2*mod(x0+elx/2,elx) - x0 - x(sel);
            dx(~sel) = x(~sel) - x0;
            %y
            sel = abs(y-y0) > ely/2;
            dy(sel) = 2*mod(y0+ely/2,ely) - y0 - y(sel);
            dy(~sel) = y(~sel) - y0;
    distE = sqrt(dx.^2+dy.^2);
    plot(pex(pick),pey(pick),'.r');
    % E<-I
    pick = m(p.nv1e+1:p.nv1,i)>0;
    x = pix(pick);
    y = piy(pick);
    npick = sum(pick);
    dx = zeros(npick,1);
    dy = zeros(npick,1);
    % enforce periodic boundary condition
            %x
            sel = abs(x-x0) > ilx/2;
            dx(sel) = 2*mod(x0+ilx/2,ilx) - x0 - x(sel);
            dx(~sel) = x(~sel) - x0;
            %y
            sel = abs(y-y0) > ily/2;
            dy(sel) = 2*mod(y0+ily/2,ily) - y0 - y(sel);
            dy(~sel) = y(~sel) - y0; 
    distI = sqrt(dx.^2+dy.^2);
    plot(x,y,'.b');
    plot(x0,y0,'*g');
    
    title([num2str(i),'th neuron, Exc, nLGN = ',num2str(nLGN(i,1))]);
    xlim([min(min(pex),min(pix)),max(max(pex),max(pix))]);
    ylim([min(min(pey),min(piy)),max(max(pey),max(piy))]);
    axis equal
%   subplot(2,3,3);
    subplot(2,2,2);
    dr = sqrt((p.dE*p.ev1x/2)^2+(p.dE*p.ev1y/2)^2);
    dC = dr/10;
    centers = (0:dC:(dr-dC)) + 0.5*dC;
    cE = hist(distE,centers);
    cI = hist(distI,centers);
    bar(centers,[cI', cE']);
    legend({'inh','exc'});
    xlabel('\mu m');
    ylabel('# of connection');
    
    % I;
    i = eg(j,2);
    subplot(2,2,3);
    hold on
%     x0 = p.iv1x;
%     y0 = p.iv1y;
%     i = (x0-1)*p.iv1y+y0;
    y0 = piy(i);
    x0 = pix(i);   
    % I<-E
    pick = m(1:p.nv1e,p.nv1e+i)>0;
    retro = mei(i,pick)>0;
    x = pex(pick);
    y = pey(pick);
    npick = sum(pick);
    dx = zeros(npick,1);
    dy = zeros(npick,1);
    % enforce periodic boundary condition
            %x
            sel = abs(x-x0) > elx/2;
            dx(sel) = 2*mod(x0+elx/2,elx) - x0 - x(sel);
            dx(~sel) = x(~sel) - x0;
            %y
            sel = abs(y-y0) > ely/2;
            dy(sel) = 2*mod(y0+ely/2,ely) - y0 - y(sel);
            dy(~sel) = y(~sel) - y0;
    distE = sqrt(dx.^2+dy.^2);
    plot(x,y,'.r');
    % I<-I
    pick = m(p.nv1e+1:p.nv1,p.nv1e+i)>0;
    x = pix(pick);
    y = piy(pick);
    npick = sum(pick);
    dx = zeros(npick,1);
    dy = zeros(npick,1);
    % enforce periodic boundary condition
            %x
            sel = abs(x-x0) > ilx/2;
            dx(sel) = 2*mod(x0+ilx/2,ilx) - x0 - x(sel);
            dx(~sel) = x(~sel) - x0;
            %y
            sel = abs(y-y0) > ily/2;
            dy(sel) = 2*mod(y0+ily/2,ily) - y0 - y(sel);
            dy(~sel) = y(~sel) - y0; 
    distI = sqrt(dx.^2+dy.^2);
    plot(x,y,'.b');
    plot(x0,y0,'*g');
    title([num2str(i),'th neuron, Inh, nLGN = ',num2str(nLGN(i,1))]);
    xlim([min(min(pex),min(pix)),max(max(pex),max(pix))]);
    ylim([min(min(pey),min(piy)),max(max(pey),max(piy))]);   
    axis equal
%   subplot(2,3,6);
    subplot(2,2,4);
    cE = hist(distE,centers);
    cI = hist(distI,centers);
    cR = hist(distE(retro),centers);
    bar(centers,[cI', cE',cR']);
    xlabel('\mu m');
    ylabel('# of connection');
    legend({'inh','exc','reciprocal'});
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hEx,['conMatFigures/eg',num2str(j),'-',matName,'.',format],printDriver,dpi);
        saveas(hEx,['conMatFigures/eg',num2str(j),'-',matName,'.fig']);
    end
end
end
%%
hPre = figure;
con.preEE = mean(sum(mee>0));
con.preEI = mean(sum(mei>0));
con.preIE = mean(sum(mie>0));
con.preII = mean(sum(mii>0));
subplot(2,2,1);
hist(sum(mee>0));
title(['pre EE: ', num2str(con.preEE)]);
xlabel('# presynaptic')
ylabel('# cell')
subplot(2,2,2);
hist(sum(mei>0));
title(['pre EI: ', num2str(con.preEI)]);
xlabel('# presynaptic')
ylabel('# cell')
subplot(2,2,3);
hist(sum(mie>0));
title(['pre IE: ', num2str(con.preIE)]);
xlabel('# presynaptic')
ylabel('# cell')
subplot(2,2,4);
hist(sum(mii>0));
title(['pre II: ', num2str(con.preII)]);
xlabel('# presynaptic')
ylabel('# cell')
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    print(hPre,['conMatFigures/pre-',matName,'.',format],printDriver,dpi);
end

disp(matName);
disp('pre:');
disp(['EE:',num2str(con.preEE),', EI:',num2str(con.preEI)]);
disp(['IE:',num2str(con.preIE),', II:',num2str(con.preII)]);

con.postEE = mean(sum(mee>0,2));
con.postEI = mean(sum(mei>0,2));
con.postIE = mean(sum(mie>0,2));
con.postII = mean(sum(mii>0,2));
hPre = figure;
subplot(2,2,1);
hist(sum(mee>0,2));
title(['post EE: ', num2str(con.postEE)]);
xlabel('# postsynaptic')
ylabel('# cell')
subplot(2,2,2);
hist(sum(mei>0,2));
title(['post EI: ', num2str(con.postEI)]);
xlabel('# postsynaptic')
ylabel('# cell')
subplot(2,2,3);
hist(sum(mie>0,2));
title(['post IE: ', num2str(con.postIE)]);
xlabel('# postsynaptic')
ylabel('# cell')
subplot(2,2,4);
hist(sum(mii>0,2));
xlabel('# postsynaptic')
ylabel('# cell')
title(['post II: ', num2str(con.postII)]);
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    print(hPre,['conMatFigures/post-',matName,'.',format],printDriver,dpi);
end
disp('post:');
disp(['EE:',num2str(con.postEE),', EI:',num2str(con.postEI)]);
disp(['IE:',num2str(con.postIE),', II:',num2str(con.postII)]);
save(['Neff-',matName,'.mat'],'con');
%%
% if reciprocal > 0.0
    hRetro = figure;
    subplot(2,2,1);
    neff0ie = sum((mie.*mei'));
    retroPer = neff0ie./sum(mie>0) * 100;
    hist(retroPer);
    xlabel('reciprocal percentage %');
    ylabel('# Cell');
    title(['mean: ',num2str(mean(retroPer)),'%']);
    subplot(2,2,2);
    hist([neff0ie' sum(mie>0)']);
    xlabel('# preIE Cell');
    legend({'reciprocal','total'});
    ylabel('# Post Cell');

    subplot(2,2,3);
    neff0ei = sum((mei.*mie'));
    retroPer = neff0ei./sum(mei>0) * 100;
    hist(retroPer);
    xlabel('reciprocal percentage %');
    ylabel('# Cell');
    title(['mean: ',num2str(mean(retroPer)),'%']);
    subplot(2,2,4);
    hist([neff0ei' sum(mei>0)']);
    xlabel('# preEI Cell');
    legend({'reciprocal','total'});
    ylabel('# Post Cell');

    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hRetro,['conMatFigures/retro-',matName,'.',format],printDriver,dpi);
    end
% end
%%
hCon = figure;
nLGNmin = min(nLGN(1:p.nv1e,1));
nLGNmax = max(nLGN(1:p.nv1e,1));
subplot(2,1,1);
hold on
f = nLGN(1:p.nv1e,1);
g = (sum(mee>0))';
hist(round(f + g/ratio),nLGNmax-nLGNmin);
plot(totalCon,0,'*r');
plot(mean(round(f+g/ratio)),0,'*g');
legend({'totalCon','prescribed','mean'});
title(['TC excitation would takes: ',num2str(100*mean(f)/(mean(g)/ratio+mean(f)),'%2.1f'),'%']);
xlabel('scaled total connections (TC+CC/ratio)');
ylabel('# Exc V1 Cells')

subplot(2,1,2);
pick = p.nv1e+(1:p.nv1i);
f = nLGN(pick,1);
g = (sum(mie>0))';
hist(round(f + g/ratio),max(nLGN(pick,1))-min(nLGN(pick,1)));
% title(['TC excitation would takes: ',num2str(100*mean(f)/(mean(g)/ratio+mean(f)),'%2.1f'),'%']);
xlabel('scaled total connections (TC+CC/ratio)');
ylabel('# Inh V1 Cells')
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    print(hCon,['conMatFigures/totalNconEff-',matName,'.',format],printDriver,dpi);
end
%load handel
%sound(y,Fs)
end
function [X, Y, Z] = multiLGNspatialKernel(posx,posy,x,y,n,s,scale)
    Aa = 14.88 * (180/pi)^2;
    Ab = Aa*0.97;
    rc = 5.61; %deg
    rs = 16.98*rc/5.61; %deg
    siga = rc/sqrt(2);
    sigb = rs/sqrt(2);
    Z = zeros(length(y),length(x));
    [X, Y] = meshgrid(x,y);
    for i = 1:n
        Z = Z + s(i)*(Aa/(siga^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/siga^2/2)...
            - Ab/(sigb^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/sigb^2/2));
    end
%     for i = 1:n
%         Z = Z + 1*(Aa/(siga^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/siga^2/2)...
%             - Ab/(sigb^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/sigb^2/2));
%     end
    Z = Z*scale;
end
