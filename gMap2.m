%% setup
drawSpecific = false;
p = struct;
p.YoudensIndex = 0.32;
p.profile = 'u'         % 'u' for uniform(single-value) profile for LGN connection strength, 'g' for gaussian
p.name = 'ndi305-40';
enorm = 0.305; 
inorm = 0.40; 
enormStd = 0.1016;
inormStd = 0.1016;
p.h = 'png';
p.esigma0 = [10.5,0.1];
p.isigma0 = [10.5,0.1];
%p.isigma0 = p.esigma0 * 1.1;
p.eAspectR0 = [1.2;0.012];
p.iAspectR0 = [1.4;0.014];
%p.iAspectR0 = p.eAspectR0*1.19;
p.se = 1;  % connection strength
%p.si = 1;
p.si = 4/3.2875;
%p.si = 1;    %l
%p.si = 1;       %b
aliasing = true;
p.pCRF = 0.0; % percentage of pure complex cells, full overlap LGN input
pMono = 0;
ep = [0, pMono,	1-p.pCRF]; % 27% 3-Stripes Niell & Stryker 2008
FontSize = 14;
set(0,'DefaultAxesFontSize',FontSize);
pPosition = [18, 180, 1200, 900];
if ~isempty(p.h)
    if strcmp(p.h,'psc2')
        printDriver = ['-de',p.h];
        p.h = 'eps';
    else
        printDriver = ['-d',p.h];
    end
    dpi = '-r150';
end
multi = 1.2;
p.lgnx = 16;%*multi;
p.lgny = 16;%*multi;
% p.lgnx = 6;
% p.lgny = 6;
% p.uniformStrength = false;
p.uniformStrength = true;
p.ev1x = 60*multi; % column 60
p.ev1y = 100*multi;	% row 100
p.iv1x = 30*multi;    % 30
p.iv1y = 50*multi;    % 50

p.dE = 6.515;
p.dI = 13.03;
% p.lgnmax = 30;  % per square degree % 0.037 per squared lgn RF, 0.029 per circled LGN RF
p.lgn_offset = 0.02; % lgn position offset
p.v1_offset = 0.02;	% v1 position offset relative to distance to nearest neuron % not implemented
p.example = 32; % examples bigger than 16 will not be displayed.
p.sigmaTheta = 0/180*pi; % orientation fluctuation of subregion
% guassian parameters, to use cutoff uniform, put seed to -1;
p.dlgn = 5/(5/180*pi)^2; % corresponding to lgns per 25 square degree
p.seed = 911; % keep it fix, data file naming system doesn't have it.
rng(p.seed);
p.average = 3;
% guassian end
p.cbounde = sqrt(2);
p.cboundi = p.cbounde*1.1;
% p.cboundi = p.cbounde*1.0;
p.sbounde = 1;
p.sboundi = 1;
%p.si = 1;
% p.boundi = sqrt(2)+0.5;
% p.bound = 2*sqrt(2*log(2));

% p.pick = true; % only use picked lgn, change jlgn in adg_INPUT for FASTER lgn input rendering.
p.pick = false;
% p.drawOrientation = true;
p.drawOrientation = false;
% mouse visual space 
%	altitude:	-50,25
%	azimuth:	-25,150
%	biggest RF:	35 degrees in diameter
%	smallest RF: 10 degrees in diameter
%	percentage of ON/OFF dominating, 0.27 for layer2/3, 0.73 for layer4
% pMono = 0.4167; % combined Layer 2 and 3
% pMono = 0.27; % Layer 2 and 3 only
%ep = [0,	pMono*1,	0.73]; % 27% 3-Stripes Niell & Stryker 2008
% p.YoudensIndex = 0.5; % 0.32 for Liu et al 2010 online method
ip = [0,    0,  1];
% p.nSubTypeE = [1,2,2,3];
% p.nSubTypeE = [1,2,2];
% p.Etypes = {' ON-OFF',' ORF',' SRF',' 3-Stripes'};
% p.Etypes = {' ON-OFF',' ORF',' SRF'};
p.Etypes = {' ORF',' SRF',' CRF'};
p.nSubTypeE = [2,2,2];
p.Etypes = strcat('Exc',p.Etypes);
p.ntypeE = length(p.Etypes);
p.Itypes = {' ORF'};
p.nSubTypeI = [2];
p.Itypes = strcat('Inh',p.Itypes);
p.ntypeI = length(p.Itypes);
p.pOn = 0.85;
core = 26.06*multi;
% core = 5;
extend = 20;
p.multi = 1;

azimuth = [-25,150]./180.*pi;
altitude = [-50,25]./180.*pi;
p.fdr = ['lgn2v1map'];
if ~exist(p.fdr,'dir')
    mkdir(p.fdr);
end

p.pinwheel = [];
p.rlgn = 0.1555;
% p.pinwheel = [2,2];
% auto data
range = (p.multi*core+2*extend)./180.*pi;
brim = (diff(altitude)-range)/2;
p.lgn_altitude = [altitude(1)+brim,altitude(2)-brim];
brim = (diff(azimuth)-range)/2;
p.lgn_azimuth = [azimuth(1)+brim,azimuth(2)-brim];

range = (p.multi*core)./180.*pi;
brim = (diff(altitude)-range)/2;
p.v1_altitude = [altitude(1)+brim,altitude(2)-brim];
brim = (diff(azimuth)-range)/2;
p.v1_azimuth = [azimuth(1)+brim,azimuth(2)-brim];

p.nlgn = p.lgnx*p.lgny;

p.nv1i = p.iv1x*p.iv1y;
p.nv1e = p.ev1x*p.ev1y;
p.nv1 = p.nv1e + p.nv1i;

%% distribute RF in exc neuron
p_subregion = rand([p.nv1e,1]);
p.eSubregion = zeros(p.nv1e,1);
p.typeE = zeros(p.nv1e,1);
mostSubregion = length(ep)-1;
 % here sigma is only FWHM, later transformed to sigma
p.esigma = zeros(p.nv1e,mostSubregion);
if mostSubregion>1
	p.ePeakDistance = zeros(p.nv1e,mostSubregion-1);
end
p.eAspectRatio = zeros(p.nv1e,mostSubregion);
p.enormDistance = zeros(p.nv1e,1);
if pMono > 0
%%%%%%%%%%%%% ON-OFF dominating
ir = (p_subregion > ep(1) & p_subregion <= ep(2));
nir = sum(ir);
p.eSubregion(ir) = 1;

% on-off axis
p.esigma(ir,1) = (10.5 + randn([nir,1])*0.1)./180.*pi;
% orth
% p.eAspectRatio(ir,1) = (1.5 + randn([nir,1])*0.1);
temp = (p.eAspectR0(1) + randn([nir,1])*p.eAspectR0(2));
temp(temp<1) = 1;
p.eAspectRatio(ir,1) = temp;
p.typeE(ir) = 1;
end
%%%%%%%%%%%%%% SRF/ORF
    ir = (p_subregion > ep(2) & p_subregion <= ep(3));
    nir = sum(ir);
    p.eSubregion(ir) = 2;
    
    % on-off axis
    p.esigma(ir,1) = (p.esigma0(1) + randn([nir,1])*p.esigma0(2))./180.*pi;
    p.esigma(ir,2) = p.esigma(ir,1);
    % peak distance
    % normDistance(ir) = 0.1 + rand(nir,1)*0.75;
    temp = zeros(nir,1)-1;
    % figure(100);
    nnn = 0;
    % nn0 = 0;
    % binranges = -0.2:0.05:1;
    % nbins = length(binranges);
    % binranges(nbins) = 1.1;
    if aliasing
        temp = enorm + randn(nir-nnn,1)*enormStd;
        neg = temp<0;
        temp(neg) = -temp(neg);
    else
        while nnn < nir
            unpicked = temp<0;
            temp(unpicked) = enorm + randn(nir-nnn,1)*enormStd;
        %     nn0 = nn0 + nir - nnn;
            temp(temp<0.0) = -1;
            nnn = nnn + sum(temp>0);
        end
    end
    % b = histc(temp,binranges);
    % b = b(1:nbins-1);
    % 
    % b0 = histc(0.3+randn(nn0,1)*0.15,binranges);
    % b0 = b0(1:nbins-1);
    % plot(binranges(1:nbins-1),[b,b0],'*');
    
    % temp(temp<0) = 0;
    
    p.enormDistance(ir) = temp; 
    % normDistance(ir) = 0.0 + rand(nir,1)*(YoudensIndex/perORF);
    p.ePeakDistance(ir,1) = 0.5*p.enormDistance(ir) .* (p.esigma(ir,1)+p.esigma(ir,2));
    
    % orth
    % p.eAspectRatio(ir,1) = (1.2 + randn([nir,1])*0.07);
    % p.eAspectRatio(ir,2) = p.eAspectRatio(ir,1);
    temp = (p.eAspectR0(1) + randn([nir,1])*p.eAspectR0(2));
    temp(temp<1) = 1.0;
    p.eAspectRatio(ir,1) = temp;
    p.eAspectRatio(ir,2) = temp;

    if pMono>0
        p.typeE(ir & p.enormDistance > p.YoudensIndex) = 3;
        p.typeE(ir & p.enormDistance <= p.YoudensIndex) = 2;
    else
        p.typeE(ir & p.enormDistance > p.YoudensIndex) = 2;
        p.typeE(ir & p.enormDistance <= p.YoudensIndex) = 1;
    end

%%%%%%%%%%%%%% CRF
if p.pCRF>0
    ir = (p_subregion > ep(3) & p_subregion <= 1);
    nir = sum(ir);
    p.eSubregion(ir) = 2;
    
    % on-off axis
    p.esigma(ir,1) = (p.esigma0(1) + randn([nir,1])*p.esigma0(2))./180.*pi;
    p.esigma(ir,2) = p.esigma(ir,1);
    % peak distance
    p.enormDistance(ir) = 0.0;
    p.ePeakDistance(ir,1) = 0.5*p.enormDistance(ir) .* (p.esigma(ir,1)+p.esigma(ir,2));
    
    % orth
    temp = (p.eAspectR0(1) + randn([nir,1])*p.eAspectR0(2));
    temp(temp<1) = 1.0;
    p.eAspectRatio(ir,1) = temp;
    p.eAspectRatio(ir,2) = temp;
    p.typeE(ir) = 3;
end

%% distribute RF in inh neuron
p_subregion = rand([p.nv1i,1]);
p.iSubregion = zeros(p.nv1i,1);
p.typeI = zeros(p.nv1i,1);
mostSubregion = length(ip)-1;
p.isigma = zeros(p.nv1i,mostSubregion);
if mostSubregion>1
	p.iPeakDistance = zeros(p.nv1i,mostSubregion-1);
end
p.iAspectRatio = zeros(p.nv1i,mostSubregion);
p.inormDistance = zeros(p.nv1i,1);

ir = (p_subregion > ip(2) & p_subregion <= ip(3));
nir = sum(ir);
p.iSubregion(ir) = 2;

% on-off axis
% p.isigma(ir,:) = (10.5 + randn([nir,2])*1)./180.*pi;
p.isigma(ir,1) = (p.isigma0(1) + randn([nir,1])*p.isigma0(2))./180.*pi;
% p.isigma(ir,1) = (12.0 + randn([nir,1])*0.1)./180.*pi;
p.isigma(ir,2) = p.isigma(ir,1);
% p.isigma(ir,2) = p.isigma(ir,1).*(1.0 + randn([nir,1])*0.07);

% orth
% p.iAspectRatio(ir,1) = (1.2 + randn([nir,1])*0.1);
% p.iAspectRatio(ir,2) = p.iAspectRatio(ir,1);
temp = (p.iAspectR0(1) + randn([nir,1])*p.iAspectR0(2));
% temp = (1.5 + randn([nir,1])*0.012);
temp(temp<1) = 1.0;
p.iAspectRatio(ir,1) = temp;
p.iAspectRatio(ir,2) = p.iAspectRatio(ir,1);

% peak distance
% temp = 0.50 + randn(nir,1)*0.1;
% temp(temp<0) = 0;
temp = zeros(nir,1)-1;
nnn = 0;
while nnn < nir
    unpicked = temp<0;
    temp(unpicked) = inorm + randn(nir-nnn,1)*inormStd;
    temp(temp<0.0) = -1;
    nnn = nnn + sum(temp>0);
end
p.inormDistance(ir) = temp;
p.iPeakDistance(ir,1) = 0.5*p.inormDistance(ir) .* (p.isigma(ir,1)+p.isigma(ir,2));

p.typeI(ir) = 1;
%%
h = figure;
subplot(2,2,1);
hist(p.eSubregion,1:length(ep)-1);
xlabel('Subregions');
ylabel('# neurons');

ORF_SRF = p.eSubregion==2;
% p.normDistance = p.ePeakDistance(ORF_SRF,1)*2 ./(p.esigma(ORF_SRF,1)+p.esigma(ORF_SRF,2));

n2sub = sum(p.eSubregion == 2);
pORF = sum(p.typeE == 1)/n2sub;
pSRF = sum(p.typeE == 2)/n2sub;
pMono = sum(p.eSubregion == 1)/p.nv1e;
pS_ORF = sum(p.eSubregion == 2)/p.nv1e;
% p3Stripes = sum(p.eSubregion == 3)/p.nv1e;
disp(['prescribed: ',num2str(pORF*100,'%3.1f'),'% ORF vs ',num2str(pSRF*100,'%3.1f'),'% SRF']);
title({['prescribed: ',num2str(pORF*100,'%3.1f'),'% ORF vs ',num2str(pSRF*100,'%3.1f'),'% SRF'], ...
        [num2str(pMono*100,'%3.1f'),'% ON/OFF, ',num2str(pS_ORF*100,'%3.1f'),'% S/ORF, ']});%,num2str(p3Stripes*100,'%3.1f'),'% 3-Stripes']});
subplot(2,2,2);
hold on
m = zeros(length(ep)-1,1);
s = zeros(length(ep)-1,1);
mp = zeros(length(ep)-2,1);
sp = zeros(length(ep)-2,1);
for k=1:(length(ep)-1)
    pick = p.eSubregion == k;
    if k > 0
        m(k) = mean(mean(p.esigma(pick,1:k).*(1+p.eAspectRatio(pick,1:k))./2,2)*180/pi);
        s(k) = std(mean(p.esigma(pick,1:k).*(1+p.eAspectRatio(pick,1:k))./2,2)*180/pi);
        if k > 1
            mp(k-1) = mean(mean(p.ePeakDistance(pick,1:k-1),2))*180/pi;
            sp(k-1) = std(mean(p.ePeakDistance(pick,1:k-1),2))*180/pi;
        end
    end
end
errorbar(1:(length(ep)-1),m,s,'^b');
errorbar((2:(length(ep)-1))+0.2,mp,sp,'*r');
legend({'RFsize','PeakDistance'},'Location','best');

subplot(2,2,3);
hist(p.iSubregion,1:length(ip)-1);
xlabel('Subregions');
ylabel('# neurons');
subplot(2,2,4);
hold on
m = zeros(length(ip)-1,1);
s = zeros(length(ip)-1,1);
mp = zeros(length(ip)-2,1);
sp = zeros(length(ip)-2,1);
for k=1:(length(ip)-1)
    pick = p.iSubregion == k;
    if k > 0
        m(k) = mean(mean(p.isigma(pick,1:k).*(1+p.iAspectRatio(pick,1:k))./2,2)*180/pi);
        s(k) = std(mean(p.isigma(pick,1:k).*(1+p.iAspectRatio(pick,1:k))./2,2)*180/pi);
        if k > 1
            mp(k-1) = mean(mean(p.iPeakDistance(pick,1:k-1),2)*180/pi);
            sp(k-1) = std(mean(p.iPeakDistance(pick,1:k-1),2)*180/pi);
        end
    end
end
errorbar(1:(length(ip)-1),m,s,'^b');
errorbar((2:(length(ip)-1))+0.2,mp,sp,'*r');
legend({'FWHM','PeakDistance'},'Location','best');
axis auto
if ~isempty(p.h)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(p.h,'fig')
        saveas(h,['lgn2v1map/',p.name,'-presets.',p.h]);
    else
        print(h,['lgn2v1map/',p.name,'-presets.',p.h],printDriver,dpi);
    end
end
%% from FWHM to sigma
p.esigma = p.esigma./(2*sqrt(2*log(2)));
p.isigma = p.isigma./(2*sqrt(2*log(2)));

[p,etheta, itheta, sp, nLGN, nSubLGN, v1Map, LGNpos, lgnStrength] = lgn2v1Map_beta(p);
save([p.name,'.mat'],'p','etheta','itheta','sp','nLGN','nSubLGN','v1Map','LGNpos', 'lgnStrength');
meanLGNe = mean(nLGN(1:p.nv1e,1));
meanLGNi = mean(nLGN(p.nv1e+(1:p.nv1i),1));
disp(['nLGN i over e = ', num2str(meanLGNi*p.si/(meanLGNe*p.se))]);
disp(p.name);
