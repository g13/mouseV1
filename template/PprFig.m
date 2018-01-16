clear all
test000 = 0;
test001 = 1;
test002 = 1;
rr = 4;
opfdr = 'PprFigb';
if rr == 5
    theme = 'ndr305-50';
    type = 'T8';
    fdr0 = 'sample-5t8021';
end
if rr == 4
    theme = 'ndi305-40';
    type = 't06';
    fdr0 = 'is40-t61l';
    fdr0f = 'isample40-tf1l';
end
lgnfile = ['1xu-',theme,'-s911'];
%coMat = ['coMa-',lgnfile];
coMat = ['coMa-60x60-',lgnfile,'_more'];
conMat = ['conMat-',theme,'_',type,'.mat'];
sProfile = ['logNormalProfile-',theme,'_',type];
%format = 'png';    
format = 'psc2';
boxOnOff = 'on';
ntheta = 12;
dtheta = pi/ntheta;
spread = false;
data0 = [fdr0,'-tcData-x4.mat'];
dataf = [fdr0f,'-tcData-x4.mat'];
test = false;
testSetup = false;
%testResult = false;
%testAnalysis = false;
%testSetup = true;
testResult = false;
testAnalysis = false;


    SetupW = 612;
    SetupH = 0.55 * SetupW;
    SetupMW0 = 0.1*SetupW;
    SetupMH0 = 0.1*SetupH;
    SetupMW1 = SetupMW0;
    SetupMH1 = SetupMH0;

    ResultW = 612;
    ResultH = 1*ResultW;
    ResultMW0 = 0.1*ResultW;
    ResultMH0 = 0.1*ResultH;
    ResultMW1 = ResultMW0;
    ResultMH1 = ResultMH0;

    AnalysisW = 612;
    AnalysisH = 0.8*AnalysisW;
    AnalysisMW0 = 0.1*AnalysisW;
    AnalysisMH0 = 0.1*AnalysisH;
    AnalysisMW1 = AnalysisMW0;
    AnalysisMH1 = AnalysisMH0;

contrastLevel=4;
conLabel = cell(contrastLevel,1);
thres = 0.0;
FWHM = 2*sqrt(2*log(2));
LGN.Aa = 15.88 * (180/pi)^2;
LGN.Ab = LGN.Aa*0.97;
rc = 5.61; %deg
rs = 16.98*rc/5.61; %deg
rc = rc * pi/180;
rs = rs * pi/180;
LGN.siga2 = rc^2/2;
LGN.sigb2 = rs^2/2;

ttt = 0:0.1:ceil(2*pi*10)/10;
rfx = @(x0,a,b,theta,t,aratio,bratio) x0 + (aratio.*a.*cos(t).*cos(theta)-bratio.*b.*sin(t).*sin(theta));
rfy = @(y0,a,b,theta,t,aratio,bratio) y0 + (aratio.*a.*cos(t).*sin(theta)+bratio.*b.*sin(t).*cos(theta));
load(lgnfile);
load(sProfile);
load(coMat,'RFcorrMat','dThetaMat');
    RFcorrMatE = RFcorrMat(1:p.nv1e,1:p.nv1e);
    RFcorrMatI = RFcorrMat(p.nv1e +(1:p.nv1i),1:p.nv1e);
    dThetaE = dThetaMat(1:p.nv1e,1:p.nv1e);
    dThetaI = dThetaMat(p.nv1e+(1:p.nv1i),1:p.nv1e);
fid = fopen(conMat);
    m = fread(fid,[p.nv1,p.nv1],'int8');
    mee = m(1:p.nv1e,1:p.nv1e);
    mei = m(p.nv1e+(1:p.nv1i),1:p.nv1e);
    mie = m(1:p.nv1e,p.nv1e+(1:p.nv1i));
    mii = m(p.nv1e+(1:p.nv1i),p.nv1e+(1:p.nv1i));
fclose(fid);
ABCLabelFontSize = 12;
tPos = [-0.1 1.0];
pPosition = [0, 0, 1280, 720];
LineWidth = 2;
set(groot,'defaultLineLineWidth',LineWidth);
set(groot,'defaultErrorbarLineWidth',LineWidth);
FontSize = 10;
txtSize = 8;
set(groot,'defaultAxesFontSize',FontSize);
set(groot,'defaultTextFontSize',txtSize);
legendSize = 8;
set(groot,'defaultLegendFontSize',legendSize);
if ~isempty(format)
    if strcmp(format,'psc2')
        printDriver = ['-de',format];
        format = 'eps';
    else
        printDriver = ['-d',format];
    end
    dpi = '-r300';
end
if testSetup
%%  Simulation Setup
%                
%   1       3        6  
%           4        7 
%   2       5        8
%
%   1. LGN: 1.1 RF_E, 1.2 RF_I, 1.3 On-off
%   2. Salt-pepper (ellipses)
%   3. Con. Prob. vs Distance 
%   4. Con. Prob. vs Orientation Diff
%   5. Con. Prob. vs RFcoeff 
%   6. E-E logNormal 
%   7. EPSP RFcoeff sample;
%   8. Con. CDF, Exc. CDF vs RFcoeff

% 1. LGN 
% 1.1   
%       1.2
% 1.3
%       1.4
%
%hLGN = figure;
hSetup = figure;
    %set(gcf,'PaperType','usletter');
    %set(gcf,'Units', 'normalized');
    %set(gcf,'Position',[0.1,0.1,0.8,0.45]);
    %set(gcf,'Position',[0.0,0.0,1.0,1.0]);

%   1.1 RF_E
LW = 4;
dxy = 24;
redge = 0.5;
%seed = 1343343;
seed = 1343145;
rng(seed);
i = randi(p.nv1e,1);
%disp(i);
ax = subplot(4,6,1);
    % x, y, width, height
    % aPos = get(ax,'Position')
    %set(ax,'OuterPosition',[0.1,0.73,0.14,0.14]);
    set(ax,'OuterPosition',[0.08,0.75,0.16,0.16]);
    hold on
    nSubregion = p.eSubregion(i);
    xlimit = [inf,-inf];
    ylimit = xlimit;
    for k = 1:nSubregion
        xtmp = [min(LGNpos(abs(v1Map{i,k}),1)) - rs*redge, max(LGNpos(abs(v1Map{i,k}),1)) + rs*redge];
        if xtmp(1) < xlimit(1), xlimit(1) = xtmp(1); end
        if xtmp(2) > xlimit(2), xlimit(2) = xtmp(2); end
        ytmp = [min(LGNpos(abs(v1Map{i,k}),2)) - rs*redge, max(LGNpos(abs(v1Map{i,k}),2)) + rs*redge];
        if ytmp(1) < ylimit(1), ylimit(1) = ytmp(1); end
        if ytmp(2) > ylimit(2), ylimit(2) = ytmp(2); end
    end
    x = xlimit(1):rs/dxy:xlimit(2);
    y = ylimit(1):rs/dxy:ylimit(2);
    for k = 1:nSubregion
        if nSubLGN{i}(k) > 0
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(v1Map{i,k},:),x,y,nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','-','LineColor','r','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','r','LineWidth',LW/2);
        else
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(-v1Map{i,k},:),x,y,-nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','-','LineColor','b','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','b','LineWidth',LW/2);
        end
    end
    axis equal
    txt = text(0,1.22,'Exc','Units','normalized','FontSize',FontSize);
    t = title('A','Units','normalized','Position',[-0.8,1.34],'FontSize',ABCLabelFontSize);
    axis(ax,'off');

%   1.2 RF_I
i = p.nv1e + randi(p.nv1i,1);
%disp(i)
ax = subplot(4,6,2);
    %aPos = get(ax,'Position')
    ip = get(ax,'Position');
    %disp(ip);
    set(ax,'Position',[ip(1)-0.03,ip(2)-0.01,ip(3),ip(4)]);
    %set(ax,'OuterPosition',[0.21,0.67,0.14,0.20]);
    hold on
    nSubregion = p.iSubregion(i-p.nv1e);
    xlimit = [inf,-inf];
    ylimit = xlimit;
    for k = 1:nSubregion
        xtmp = [min(LGNpos(abs(v1Map{i,k}),1)) - rs*redge, max(LGNpos(abs(v1Map{i,k}),1)) + rs*redge];
        if xtmp(1) < xlimit(1), xlimit(1) = xtmp(1); end
        if xtmp(2) > xlimit(2), xlimit(2) = xtmp(2); end
        ytmp = [min(LGNpos(abs(v1Map{i,k}),2)) - rs*redge, max(LGNpos(abs(v1Map{i,k}),2)) + rs*redge];
        if ytmp(1) < ylimit(1), ylimit(1) = ytmp(1); end
        if ytmp(2) > ylimit(2), ylimit(2) = ytmp(2); end
    end
    x = xlimit(1):rs/dxy:xlimit(2);
    y = ylimit(1):rs/dxy:ylimit(2);
    for k = 1:nSubregion
        if nSubLGN{i}(k) > 0
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(v1Map{i,k},:),x,y,nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','-','LineColor','r','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','r','LineWidth',LW/2);
        else
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(-v1Map{i,k},:),x,y,-nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','-','LineColor','b','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','b','LineWidth',LW/2);
        end
    end
    axis equal
    txt = text(0,1.15,'Inh','Units','normalized','FontSize',FontSize);
    axis(ax,'off');

%   1.3 Hist Norm. Distance
ax = subplot(4,6,7);
    %aPos = get(ax,'Position')
    set(ax,'OuterPosition',[0.12,0.58,0.13,0.14]);
    %md = max(p.enormDistance);
    md = 0.8;
    gray = 0.3;
    dark = 0.3;
    dbin = 0.05;
    xlimit = ceil(md/dbin)*dbin;
    edges = 0:dbin:xlimit;
    [N,~] = histcounts(p.enormDistance, edges);
    bar(edges(1:length(edges)-1)+dbin/2,N,'FaceColor',dark*ones(1,3),'EdgeColor',gray*ones(1,3));
    set(ax,'XTick',[0,0.4,0.8],'XTickLabel',{'0','0.4','0.8'},'FontSize',FontSize);
    set(ax,'YColor','None','YTick',[]);
    %ylabel('# Exc');
    xlabel('norm. distance','FontSize',FontSize);
    xlim([0,xlimit]);
    box(ax,'off');

%   1.4 On-Off.png
ax = subplot(4,6,8);
    img = imread('On-Off.png');
    image(img);
    axis(ax,'off');
    %disp(get(ax,'Position'));
    set(ax,'Position',[0.27,0.62,0.060,0.102]);
    %disp(get(ax,'Position'));

%   2. Salt-pepper (ellipses)
    mks = 12;
    dev1y = (p.v1_altitude(2)-p.v1_altitude(1))/p.ev1y; % for correction upon periodical boundary
    axis_ev1y = linspace(p.v1_altitude(1)+0.5*dev1y,p.v1_altitude(2)-0.5*dev1y,p.ev1y);
 
    dev1x = (p.v1_azimuth(2)-p.v1_azimuth(1))/p.ev1x;
    axis_ev1x = linspace(p.v1_azimuth(1)+0.5*dev1x,p.v1_azimuth(2)-0.5*dev1x,p.ev1x);
    
    ev1posx = ones(p.ev1y,1) * axis_ev1x;
    ev1posy = axis_ev1y' * ones(1,p.ev1x);
    
    div1y = (p.v1_altitude(2)-p.v1_altitude(1))/p.iv1y;
    axis_iv1y = linspace(p.v1_altitude(1)+0.5*div1y,p.v1_altitude(2)-0.5*div1y,p.iv1y);
    
    div1x = (p.v1_azimuth(2)-p.v1_azimuth(1))/p.iv1x;
    axis_iv1x = linspace(p.v1_azimuth(1)+0.5*div1x,p.v1_azimuth(2)-0.5*div1x,p.iv1x);
    
    iv1posx = ones(p.iv1y,1) * axis_iv1x;
    iv1posy = axis_iv1y' * ones(1,p.iv1x);
    ev1pos = [reshape(ev1posx,[p.nv1e,1]),reshape(ev1posy,[p.nv1e,1])];
    iv1pos = [reshape(iv1posx,[p.nv1i,1]),reshape(iv1posy,[p.nv1i,1])];
    v1pos = [ev1pos; iv1pos];
ax = subplot(2,3,4);
    %aPos = get(ax,'Position')
    set(ax,'Position',[0.1,0.1,0.25,0.42]);
    hold on
    plot( LGNpos(:,1)*180/pi, LGNpos(:,2)*180/pi,'.k','MarkerSize',mks-2);
    hsv = zeros(p.nv1,3);
    hsv(:,1) = [etheta;itheta]./pi;
    hsv(:,2) = 0.9;
    hsv(:,3) = 0.9;
    color = hsv2rgb(hsv);
    v1pos = v1pos./pi*180;
    ndraw = ceil(0.2*p.nv1);
    drawlist = randi(p.nv1,[ndraw,1]);
    for i=1:ndraw
        j = drawlist(i);
        plot(v1pos(j,1), v1pos(j,2) ,'.','Color',color(j,:),'MarkerSize',mks); 
    end
    %example = p.nv1 - 1*p.iv1x - 2;
    example = p.nv1-4; 
    iexample = example-p.nv1e;
    %plot(v1pos(example,1), v1pos(example,2) ,'*','Color',color(i,:),'MarkerSize',mks-2); 
    plot(v1pos(example,1), v1pos(example,2) ,'*k','MarkerSize',mks-2); 
    sigma = p.isigma(iexample);
    sigmb = sigma * p.iAspectRatio(iexample,1);
    nSubregion = p.iSubregion(iexample);
    gtheta = 3*pi/4;
    %gtheta = itheta(iexample);
    for k = 1:nSubregion 
        if nSubLGN{example}(k) > 0
            peak = mean(LGNpos(v1Map{example,k},:));
            xxx2 = rfx(peak(1),sigma,sigmb,gtheta,ttt,FWHM,FWHM)*180/pi;         
            yyy2 = rfy(peak(2),sigma,sigmb,gtheta,ttt,FWHM,FWHM)*180/pi;
            plot(xxx2,yyy2,'-r','LineWidth',3); 
        else
            peak = mean(LGNpos(-v1Map{example,k},:));
            xxx2 = rfx(peak(1),sigma,sigmb,gtheta,ttt,FWHM,FWHM)*180/pi;         
            yyy2 = rfy(peak(2),sigma,sigmb,gtheta,ttt,FWHM,FWHM)*180/pi;
            plot(xxx2,yyy2,'-b','LineWidth',3); 
        end
    end
    xlim([min(LGNpos(:,1)),max(LGNpos(:,1))]*180/pi);
    ylim([min(LGNpos(:,2)),max(LGNpos(:,2))]*180/pi);
    axis equal tight

    xlabel('degree'); ylabel('degree');

    t = title('B','Units','normalized','Position',[-0.1,1.05],'FontSize',ABCLabelFontSize);

%   3. Con. Prob. vs Distance 
elx = (p.ev1x) * p.dE;
ely = (p.ev1y) * p.dE;
ilx = (p.iv1x) * p.dI;
ily = (p.iv1y) * p.dI;
pex = ones(p.ev1y,1)*((1:p.ev1x)-0.5).*p.dE;
pex = reshape(pex,[p.nv1e,1]); % id goes along y first (column major)
pey = ((1:p.ev1y)-0.5)'*ones(1,p.ev1x).*p.dE;
pey = reshape(pey,[p.nv1e,1]);

pix = ones(p.iv1y,1) * (((1:p.iv1x)-0.5).*p.dI);
pix = reshape(pix,[p.nv1i,1]);
piy = (((1:p.iv1y)-0.5)'.*p.dI) * ones(1,p.iv1x);
piy = reshape(piy,[p.nv1i,1]);

    dbin = 25;
    binranges = 0:dbin:300;
    nbins = length(binranges) - 1;
    xxxx = binranges(1:nbins)+dbin/2;

ax = subplot(3,3,2);
    hold on
    pickedNeighbor = mee>0;
    % pre E
    src.x = pex;
    src.y = pey;
    src.n = p.nv1e;
    src.lx = elx;
    src.ly = ely;
    % post E
    tar.x = pex;
    tar.y = pey;
    tar.n = p.nv1e;
    hlx = src.lx/2;
    hly = src.ly/2;
    d2 = zeros(p.nv1e,p.nv1e);
    dx = zeros(src.n,1);
    dy = zeros(src.n,1);
    for i=1:tar.n
        %x
        pick = abs(src.x-tar.x(i)) > hlx;
        dx(pick) = 2*mod(tar.x(i)+hlx,src.lx) - tar.x(i) - src.x(pick);
        %assert(sum(dx(pick) >= hlx + 1e-5)==0);
        dx(~pick) = src.x(~pick) - tar.x(i);
        %y
        pick = abs(src.y-tar.y(i)) > hly;
        dy(pick) = 2*mod(tar.y(i)+hly,src.ly) - tar.y(i) - src.y(pick);
        %assert(sum(dy(pick) >= hly + 1e-5)==0);
        dy(~pick) = src.y(~pick) - tar.y(i);

        d2(:,i) = sqrt(dx.^2 + dy.^2);
    end
    d2(~pickedNeighbor) = -1;
    [nr, ~] = histc(d2,binranges);
    nr = nr(1:nbins,:);
    nr = nr./(ones(nbins,1)*sum(nr))*100;
    mnr = mean(nr,2);
    stdnr = std(nr,1,2);
    errorbar(xxxx,mnr,stdnr,'-r');
     
    pickedNeighbor = mei>0;
    % pre I
    src.x = pix;
    src.y = piy;
    src.n = p.nv1i;
    src.lx = ilx;
    src.ly = ily;
    % post E
    tar.x = pex;
    tar.y = pey;
    tar.n = p.nv1e;
    hlx = src.lx/2;
    hly = src.ly/2;
    d2 = zeros(p.nv1i,p.nv1e);
    dx = zeros(src.n,1);
    dy = zeros(src.n,1);
    for i=1:tar.n
        %x
        pick = abs(src.x-tar.x(i)) > hlx;
        dx(pick) = 2*mod(tar.x(i)+hlx,src.lx) - tar.x(i) - src.x(pick);
        %assert(sum(dx(pick) >= hlx + 1e-5)==0);
        dx(~pick) = src.x(~pick) - tar.x(i);
        %y
        pick = abs(src.y-tar.y(i)) > hly;
        dy(pick) = 2*mod(tar.y(i)+hly,src.ly) - tar.y(i) - src.y(pick);
        %assert(sum(dy(pick) >= hly + 1e-5)==0);
        dy(~pick) = src.y(~pick) - tar.y(i);

        d2(:,i) = sqrt(dx.^2 + dy.^2);
    end
    d2(~pickedNeighbor) = -1;
    [nr, ~] = histc(d2,binranges);
    nr = nr(1:nbins,:);
    nr = nr./(ones(nbins,1)*sum(nr))*100;
    mnr = mean(nr,2);
    stdnr = std(nr,1,2);
    errorbar(xxxx,mnr,stdnr,'-b');

    xlabel('Distance \mu m');
    ylabel('% Presynaptic');
    ylim([0,inf]);
    xlim([0,binranges(end)]);
    t = title('C','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);

%   4. Con. Prob. vs Orientation Diff

ax = subplot(3,3,5);
    hold on
    pickedNeighbor = mee>0;
    binranges = 0.0:dtheta:(pi/2);
    nbins = length(binranges)-1;
    binranges(nbins+1) = pi/2+0.1;
    sel = true(p.nv1e,1);
    pickedNeighborT = pickedNeighbor(:,sel);
    dThetaMatC = dThetaE(:,sel);
    dThetaMatC(~pickedNeighborT) = binranges(1)-1;
    [binCounts,ind] = histc(dThetaMatC,binranges);
    binCounts = binCounts(1:nbins,:);
    binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
    mbinCounts = mean(binCounts,2);
    mbinCountsE = mbinCounts./sum(mbinCounts)*100;
    stdbinCountsE = std(binCounts,1,2);

    xxxx = (binranges(1:nbins)+dtheta/2)*180/pi;
    errorbar(xxxx,mbinCountsE,stdbinCountsE,'-r');

    pickedNeighbor = mei>0;
    dThetaMatC = dThetaI;
    dThetaMatC(~pickedNeighbor) = binranges(1)-1;
    binCounts = histc(dThetaMatC,binranges);
    binCounts = binCounts(1:nbins,:);
    binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
    mbinCounts = mean(binCounts,2);
    mbinCounts = mbinCounts./sum(mbinCounts)*100;
    stdbinCounts = std(binCounts,1,2);

    errorbar(xxxx,mbinCounts,stdbinCounts,'-b');

    %dThetaIE = itheta*ones(1,p.nv1e)-ones(p.nv1i,1)*etheta';
    %dThetaIE = abs(dThetaIE); 
    %pick = dThetaIE > pi/2;
    %dThetaIE(pick) = pi - dThetaIE(pick);
    %dThetaMatC = dThetaIE';

    %pickedNeighbor = mie>0;
    %dThetaMatC(~pickedNeighbor) = binranges(1)-1;
    %binCounts = histc(dThetaMatC,binranges);
    %binCounts = binCounts(1:nbins,:);
    %binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
    %mbinCounts = mean(binCounts,2);
    %mbinCounts = mbinCounts./sum(mbinCounts)*100;
    %stdbinCounts = std(binCounts,1,2);

    %errorbar(xxxx,mbinCounts,stdbinCounts,':r');

    %dThetaII = itheta*ones(1,p.nv1i)-ones(p.nv1i,1)*itheta';
    %dThetaII = abs(dThetaII); 
    %pick = dThetaII > pi/2;
    %dThetaII(pick) = pi - dThetaII(pick);
    %dThetaMatC = dThetaII';

    %pickedNeighbor = mii>0;
    %dThetaMatC(~pickedNeighbor) = binranges(1)-1;
    %binCounts = histc(dThetaMatC,binranges);
    %binCounts = binCounts(1:nbins,:);
    %binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
    %mbinCounts = mean(binCounts,2);
    %mbinCounts = mbinCounts./sum(mbinCounts)*100;
    %stdbinCounts = std(binCounts,1,2);

    %errorbar(xxxx,mbinCounts,stdbinCounts,':b');

    ylabel('% \rightarrow Exc','Interpreter','tex');
    ylim([0,inf]);
    xlim([0,90]);
    xlabel('\Delta\theta');
    t = title('D','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);

%   5. Con. Prob. vs RFcoeff 
subplot(3,3,8);
    hold on
    pickedNeighbor = mee>0;
    dbin = 0.1;
    binranges = -1:dbin:1;
    nbins = length(binranges)-1;
    binranges(nbins+1) = 1.1;
    sel = true(p.nv1e,1);

    pickedNeighborT = pickedNeighbor(:,sel);
    RFcorr = RFcorrMatE(:,sel);
    RFcorr(~pickedNeighborT) = binranges(1)-1;
    [binCounts,ind] = histc(RFcorr,binranges);
    binCounts = binCounts(1:nbins,:);
    binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
    mbinCounts = mean(binCounts,2);
    mbinCountsE = mbinCounts./sum(mbinCounts)*100;
    stdbinCountsE = std(binCounts,1,2);

    xxxx = (binranges(1:nbins)+dbin/2);
    errorbar(xxxx,mbinCountsE,stdbinCountsE,'-r');

    pickedNeighbor = mei>0;
    RFcorr = RFcorrMatI;
    RFcorr(~pickedNeighbor) = binranges(1)-1;
    binCounts = histc(RFcorr,binranges);
    binCounts = binCounts(1:nbins,:);
    binCounts = binCounts./(ones(nbins,1)*sum(binCounts,1))*100;
    mbinCounts = mean(binCounts,2);
    mbinCounts = mbinCounts./sum(mbinCounts)*100;
    stdbinCounts = std(binCounts,1,2);

    errorbar(xxxx,mbinCounts,stdbinCounts,'-b');

    %ylabel('% Exc');
    ylabel('% \rightarrow Exc','Interpreter','tex');
    ylim([0,inf]);
    xlim([-1,1]);
    xlabel('RF Correlation.');
    t = title('E','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);

%   6. E-E logNormal 
    ax = subplot(3,3,3);
        mu = 0.45;
        mu0 = 0.45;
        sigma = 1.16;
        n = 400;
        nbins = 30;
        nsig = 3;
        m = log(mu/sqrt(1+sigma^2/mu^2));
        sig = sqrt(log(1+sigma^2/mu^2));
        strengthID = zeros(n,1);
        A = 2/erfc(-(log(nsig*sig)-m)/(sig*sqrt(2)));
        logNormalCDF = @(x,m,sig) 1/2*erfc(-(log(x)-m)/(sig*sqrt(2)))*A;
        x0max = min(log(4.0),m+nsig*sig);
        x0min = max(log(5e-3),m-nsig*sig);
        x0 = exp(linspace(x0min,x0max,nbins));
        slist = x0;
        cuts = logNormalCDF(x0,m,sig);
        cuts = cuts/max(cuts);
        csnpick = round(cuts*n);
        npick = diff([0,csnpick]);
        last = 0;
        for i = 1:nbins
            if npick(i) > 0
                if last+npick(i) > n
                   npick(i) = n - last; 
                end
                ipick = last+(1:npick(i));
                strengthID(ipick) = i;
                last = last + npick(i);
            end
            if last == n
                break;
            end
        end
        
        slist = (mu0/mean(slist(strengthID)))*slist;
        logNormal = @(x,m,sig) 1./(x*sqrt(2*pi)*sig).*exp(-(log(x)-m).^2/(2*sig^2));
        pick = npick>0;
        x = slist(pick);
        y = npick(pick);
        nx = sum(pick);
        dx = zeros(1,nx);
        for i = 1:nx
            switch i
                case 1
                    dx(i) = x(i+1)-x(i);
                case nx
                    dx(i) = x(i) - x(i-1);
                otherwise
                    dx(i) = 0.5*(x(i+1)-x(i-1));
            end
            
        end
        yprob = logNormal(x,m,sig).*dx;
        yprob = yprob./sum(yprob);
        
        hold on
        ss = slist(strengthID);
        means = mean(ss);
        stds = std(ss);
        dedge = 0.2;
        edges = 0:dedge:5;
        nedge = length(edges)-1;
        gray = 0.3;
        %histogram(ss,edges,'FaceColor',gray*ones(1,3),'EdgeColor',gray*ones(1,3));
        [N,~] = histcounts(ss,edges);
        bar(edges(1:nedge)+dedge/2,N,'FaceColor',gray*ones(1,3),'EdgeColor',gray*ones(1,3));
        xlim([0,edges(end)]);
        ylabel('# Exc');
        xlabel('EPSP mV');
        t = title('F','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);

        %ax2 = axes('Position',[0.5,0.5,0.3,0.3],'Units','normalized','Box','off','XScale','log');
        pos = get(ax,'Position');
        ax2 = axes('Position',[pos(1)+0.5*pos(3),pos(2)+0.5*pos(4),pos(3)*0.5,pos(4)*0.5],'Units','normalized','Box','off','XScale','log');
        hold(ax2,'on');
        
        semilogx(x,y,'-o','LineWidth',1);
        semilogx(x+dx,yprob*n,'-','LineWidth',1);
        set(ax2,'XTick',[1e-2,1e-1,1e0,1e1],'XTicklabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}'},'XMinorTick','on','FontSize',FontSize*0.5);
        title(ax2,{'logNormal',[num2str(means,'%0.2f'),'\pm',num2str(stds,'%0.2f'),'mV']},'FontSize',FontSize*0.75);
        xlim([5e-3,10]);
%         ylim([0,inf]);
        %legend({'analytical','data'});
        
%   7. EPSP RFcoeff sample;
ax1 = subplot(3,3,6);
    p0 = get(ax1,'Position');
    %disp(p0);
    p1 = get(ax1,'OuterPosition');
    %disp(p1);
    set(ax1,'OuterPosition',[0.64,0.43,0.2869,0.2679]);
    set(ax1,'Position',[0.6916,0.45,0.2134,0.1934]);
    p0 = get(ax1,'Position');
    %disp(p0);
    slist = [1,slist];

    rng(seed+1);
    i = randi(p.nv1e,1);

    dbin =0.2;
    pickedNeighbor = mee(:,i) > 0;
    pickedNeighbor(i) = false;
    binranges = -1.0:dbin:1;
    nbins = length(binranges)-1;
    binranges(nbins+1) = 1.1;
    binCounts = histc(RFcorrMatE(pickedNeighbor,i),binranges);
    binCounts = binCounts(1:nbins);
    if spread
        %EPSPs = profiles(mee(pickedNeighbor,i),i);
        EPSPs = slist(mee(pickedNeighbor,i),i);
    else
        %EPSPs = profiles(mee(pickedNeighbor,i));
        EPSPs = slist(mee(pickedNeighbor,i));
    end
    x = binranges(1:nbins)+dbin/2;

    pos = get(ax1,'Position');
    ax2 = axes('Position',pos);
    %set(ax1,'Box','off');

    axes(ax1);
    gray = 0.5;
    b = bar(x,100*binCounts./sum(binCounts),'FaceColor',gray*ones(1,3),'EdgeColor',gray*ones(1,3)); 
    set(ax1,'Box','off','Color','none');
    
    axes(ax2);
    hold(ax2,'on');
    plot(RFcorrMatE(pickedNeighbor,i),EPSPs,'.k','MarkerSize',mks);
    set(ax2,'YAxisLocation','Right','XTick',[],'Box','off','Color','none');
    %set(ax2,'YTick',[1,round(y),100],'YTickLabel',{'1',num2str(round(y)),'100'});
    %txt = text(1.1,y,num2str(round(y)),'Units','data');

    xlim(ax1,[-1,1]);
    ylim(ax1,[0,30]);
    linkaxes([ax1,ax2],'x');
    xt = xlabel('RF correlation');
    xt.Units='normalized';
    xt.Position=[0.5,-0.17,0];
    %set(xt,'Position',get(xt,'Position')+[0,0,0.2]);
    ylabel(ax1,'% Connection');
    ylabel(ax2,'EPSP mV');
    t = title('G','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);
%   8. Con. CDF, Exc. CDF vs RFcoeff
ax1 = subplot(3,3,9);
    %disp('H:');
    p0 = get(ax1,'OuterPosition');
    %disp(p0);

    dbin0 = 0.2;
    binranges0 = -1.0:dbin0:1;
    nbins0 = length(binranges0)-1;
    [binCounts0,~] = histcounts(RFcorrMatE,binranges0);       
    binCounts0 = binCounts0./sum(binCounts0)*100;
    binranges0 = binranges0(1:nbins0) + dbin0/2;

    dbin = 0.05;
    pickedNeighbor = mee>0;
    binranges = -1.0:dbin:1;
    nbins = length(binranges)-1;
    [~,~,ind] = histcounts(RFcorrMatE,binranges);       
    ind = reshape(ind,[p.nv1e,p.nv1e]);

    RFcorrMatE(~pickedNeighbor) = binranges(1)-1;
    [binCounts,~] = histcounts(RFcorrMatE,binranges);
    binCounts = binCounts./sum(binCounts)*100;

    EPSPslice = zeros(nbins,1);
    for j = 1:nbins
        if spread
            for i = 1:p.nv1e
                ineighbor = pickedNeighbor(:,i);
                EPSPslice(j) = EPSPslice(j) + sum(profiles(mee((ind(:,i)==j)&ineighbor,i),i));
            end
            EPSPslice(j) = EPSPslice(j)/p.nv1e;
        else
            EPSPslice(j) = sum(profiles(mee((ind==j)&pickedNeighbor)))/p.nv1e;
        end
    end
    EPSPslice = EPSPslice/sum(EPSPslice)*100;
    binranges = binranges(1:nbins)+dbin/2;
    pos = get(ax1,'Position');
    ax2 = axes('Position',pos);
    %set(ax1,'Box','off');

    axes(ax1);
    gray = 0.5;
    b = bar(binranges0,binCounts0,'FaceColor',gray*ones(1,3),'EdgeColor',gray*ones(1,3)); 
    set(ax1,'Box','off','Color','none');

    axes(ax2);
    hold(ax2,'on');
    tEPSPslice = cumsum(EPSPslice);
    tbinCounts = cumsum(binCounts');
    plot(binranges,tbinCounts,'-','LineWidth',4);
    plot(binranges,tEPSPslice,'-','LineWidth',4);

    p50 = find(tEPSPslice>50,1,'first'); 
    r = (50-tEPSPslice(p50-1))/(tEPSPslice(p50) - tEPSPslice(p50-1));
    x = binranges(p50-1) + dbin * r;
    y = tbinCounts(p50-1) + r * (tbinCounts(p50)-tbinCounts(p50-1));
    
    plot([x,1],[50,50],':k','LineWidth',2);
    plot([x,x],[50,y],':k','LineWidth',2);
    plot([x,1],[y,y],':k','LineWidth',2);
    set(ax2,'YAxisLocation','Right','XTick',[],'Box','off','Color','none');
    set(ax2,'YTick',[1,50,round(y)],'YTickLabel',{'1','50',num2str(round(y))});
    %set(ax2,'YTick',[1,round(y),100],'YTickLabel',{'1',num2str(round(y)),'100'});
    %txt = text(1.1,y,num2str(round(y)),'Units','data');

    xlim(ax1,[-1,1]);
    linkaxes([ax1,ax2],'x');
    ylabel(ax1,'% PDF');
    ylabel(ax2,'% CDF');
    xlabel(ax1,'RF Correlation');
    [leg,icons,~,~] = legend(ax2,{'Prob.CDF','Cort.Exc'},'Location','NorthWest','FontSize',legendSize);
    leg.Box = 'off';
    lp = leg.Position;
    leg.Position = [lp(1)-0.023,lp(2)+0.020,lp(3),lp(4)];
    mp = icons(3).XData;
    icons(3).XData = [mp(1),mp(1)+(mp(2)-mp(1))*0.2];
    tp = icons(1).Position;
    icons(1).Position = [tp(1)-(mp(2)-mp(1))*0.8,tp(2),tp(3)];

    mp = icons(5).XData;
    icons(5).XData = [mp(1),mp(1)+(mp(2)-mp(1))*0.2];
    tp = icons(2).Position;
    icons(2).Position = [tp(1)-(mp(2)-mp(1))*0.8,tp(2),tp(3)];
    ylim(ax2,[0,100]);
    t = title('H','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);
    %txt = text(0.05,1.1,{['50% of cort. excitation by'],['top ',num2str(100-round(y)),'% correlated connection']},'Units','normalized','FontSize',txtSize,'FontWeight','bold');

saveas(hSetup,[opfdr,'/Setup-',theme,'.fig']);
if ~isempty(format)
    set(gcf,'PaperUnits','points');
    set(gcf,'PaperSize',[SetupW,SetupH]);
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperPosition',[0,0,SetupW,SetupH]);
    set(gcf,'Units','points');
    set(gcf,'OuterPosition',[SetupMW0,SetupMH0,SetupW-SetupMW1-SetupMW0,SetupH-SetupMH1-SetupMH0]);
    print(hSetup,[opfdr,'/Setup-',theme,'.',format],printDriver,dpi);
end
end
load(data0);
if test001
    figure
    rasterTC = zeros(ntheta+1,p.nv1,contrastLevel)-1;
    EmeanTC = zeros(ntheta+1,contrastLevel)-1;
    EquaLTC = zeros(ntheta+1,contrastLevel)-1;
    EquaUTC = zeros(ntheta+1,contrastLevel)-1;
    EmediTC = zeros(ntheta+1,contrastLevel)-1;
    ImeanTC = zeros(ntheta+1,contrastLevel)-1;
    IquaLTC = zeros(ntheta+1,contrastLevel)-1;
    IquaUTC = zeros(ntheta+1,contrastLevel)-1;
    ImediTC = zeros(ntheta+1,contrastLevel)-1;
    TCgLGN = rasterTC;
    TCgLGN_F1 = rasterTC;
    TCgE = rasterTC;
    TCgI = rasterTC;
    EmeanTCgLGN = EmeanTC;
    ImeanTCgLGN = ImeanTC;
    EmeanTCgLGN_F1 = EmeanTC;
    ImeanTCgLGN_F1 = ImeanTC;
    EmeanTCgE = EmeanTC;
    ImeanTCgE = ImeanTC;
    EmeanTCgI = EmeanTC;
    ImeanTCgI = ImeanTC;
    %4 norm E
    ax = subplot(1,2,1);
        ip = get(ax,'Position');
        set(ax,'Position',[ip(1)+0.03,ip(2)+0.02,ip(3)-0.09,ip(4)]);
        hold on
        maxEmeanTC = ones(ntheta+1,1) * max(EmeanTC);
        maxEmeanTCgLGN = ones(ntheta+1,1) * max(EmeanTCgLGN);
        maxEmeanTCgLGN_F1 = ones(ntheta+1,1) * max(EmeanTCgLGN_F1);
        maxEmeanTCgE = ones(ntheta+1,1) * max(EmeanTCgE);
        maxEmeanTCgI = ones(ntheta+1,1) * max(EmeanTCgI);
        for i = 1:contrastLevel
            %plot(relTheta, EmeanTCgLGN(:,i)./maxEmeanTCgLGN(:,i),LineScheme{i},'Color','g');
            plot(relTheta, EmeanTCgLGN_F1(:,i)./maxEmeanTCgLGN(:,i),LineScheme{i},'Color','m');
            plot(relTheta, EmeanTCgE(:,i)./maxEmeanTCgE(:,i),LineScheme{i},'Color','r');
            plot(relTheta, EmeanTC(:,i)./maxEmeanTC(:,i),LineScheme{i},'Color','k');
        end
        ylabel('normalized');
        xlabel('deg');
        yy = ylim();
        ylim([0,yy(2)]);
        xticks = -90:45:90;
        xlim([relTheta(1),relTheta(end)]);
        set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'));
        t = title('D','Units','normalized','Position',[-0.3,1.05],'FontSize',ABCLabelFontSize);
    %6 norm I
    ax = subplot(1,2,2);
        ip = get(ax,'Position');
        set(ax,'Position',[ip(1)+0.03,ip(2)+0.02,ip(3)-0.09,ip(4)]);
        hold on
        maxImeanTC = ones(ntheta+1,1) * max(ImeanTC);
        maxImeanTCgLGN = ones(ntheta+1,1) * max(ImeanTCgLGN);
        maxImeanTCgLGN_F1 = ones(ntheta+1,1) * max(ImeanTCgLGN_F1);
        maxImeanTCgE = ones(ntheta+1,1) * max(ImeanTCgE);
        maxImeanTCgI = ones(ntheta+1,1) * max(ImeanTCgI);
        for i = 1:contrastLevel
            %plot(relTheta, ImeanTCgLGN(:,i)./maxImeanTCgLGN(:,i),LineScheme{i},'Color','g');
            plot(relTheta, ImeanTCgLGN_F1(:,i)./maxImeanTCgLGN(:,i),LineScheme{i},'Color','m');
            plot(relTheta, ImeanTCgE(:,i)./maxImeanTCgE(:,i),LineScheme{i},'Color','r');
            plot(relTheta, ImeanTC(:,i)./maxImeanTC(:,i),LineScheme{i},'Color','k');
        end
        ylabel('norm. max');
        xlabel('deg');
        yy = ylim();
        ylim([0,yy(2)]);
        title('Inh');
        xticks = -90:45:90;
        xlim([relTheta(1),relTheta(end)]);
        set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'));
        t = title('F','Units','normalized','Position',[-0.3,1.05],'FontSize',ABCLabelFontSize);
    saveas(gcf,'test001.fig');
end
if test002
    figure;
    pos = {[1,2,7,8,13,14,19,20],[25,26,31,32]};
    columns = 2;
    nl(1) = length(pos{1});
    nl(2) = length(pos{2});
    rows = nl./columns;
    for k = 1:2
        if k==1
            neuronlist = randi(p.nv1e,[nl(k),1]);
            level = 1;
            type0 = [true(p.nv1e,1);false(p.nv1i,1)];
            types = type0(PKsortedID(:,level));
            candidates = find(types);
            ncand = length(candidates);
            neuronlist(1) = PKsortedID(candidates(1),level);
            neuronlist(3) = PKsortedID(candidates(floor(ncand/2)),level);
            neuronlist(5) = PKsortedID(candidates(ncand),level);
            level = 1;
            efired = nP(level).pkrate(dCVsortedID) > nP(level).br(dCVsortedID) & nP(level).ei(dCVsortedID) > 0.5 &...
                     nP(level).pkrate(dCVsortedID) > thres;
            type0 = [true(p.nv1e,1); false(p.nv1i,1)];
            types = type0(dCVsortedID);
            candidates = find(efired(dCVsortedID) & types); 
            ncand = length(candidates); 
            neuronlist(2) = dCVsortedID(candidates(2));
            neuronlist(4) = dCVsortedID(candidates(floor(ncand/2)));
            neuronlist(6) = dCVsortedID(candidates(ncand-1));
    
            level = 1;
            type0 = [true(p.nv1e,1);false(p.nv1i,1)];
            types = type0(PKsortedID(:,level));
            candidates = find(types);
            ncand = length(candidates);
            neuronlist(7) = PKsortedID(candidates(ncand),level);
            neuronlist(8) = PKsortedID(candidates(1),level);
            disp(neuronlist);
        else
            level = 2;
            ifired = nP(level).pkrate(dCVsortedID) > nP(level).br(dCVsortedID) & nP(level).ei(dCVsortedID) < 0.5 &...
                     nP(level).pkrate(dCVsortedID) > thres;
            neuronlist = p.nv1e + randi(p.nv1i,[nl(k),1]);
            type0 = [false(p.nv1e,1); true(p.nv1i,1)];
            types = type0(dCVsortedID);
            candidates = find(ifired(dCVsortedID) & types); 
            ncand = length(candidates); 
            neuronlist(1) = dCVsortedID(candidates(floor(ncand/1.2)));
            neuronlist(2) = dCVsortedID(candidates(floor(ncand/2)));
            neuronlist(3) = dCVsortedID(candidates(ncand));
            disp(neuronlist);
        end
        for j = 1:nl(k)
            i = neuronlist(j); 
            ax = subplot(sum(rows),columns*3,pos{k}(j));
            hold on
            inputAngle = nP(contrastLevel).priA(i);
            ipA = round(inputAngle*180/pi/dtheta)+1;
            if ipA < ntheta/2+1
                half = [(ipA+ntheta/2):ntheta,1:(ipA+ntheta/2)];
            else 
                half = [(ipA-ntheta/2):ntheta,1:(ipA-ntheta/2)];
            end
            if inputAngle >= pi/2
                pA = inputAngle - pi/2;
            else
                pA = inputAngle + pi/2;
            end
            grating_ipA = round(pA*180/pi/dtheta)+1;
            thetas = ((grating_ipA-ntheta/2-1):(grating_ipA+ntheta/2-1))*dtheta;
            for ci=1:contrastLevel
        
                rho = [tC(ci).frate(i,:),tC(ci).frate(i,1)];
                plot(thetas,rho(half),'Color','k','LineStyle',LineScheme{ci});
            end
            xlim([thetas(1),thetas(end)]);
            ylim([0,1.1*max(rho)]);
            set(ax,'XTick',[thetas(1),thetas(1)+90,thetas(end)]);
         
            if j == 1 
                if k == 1
                    ylabel('Hz');
                    xlabel('deg');
                    t = title('A','Units','normalized','Position',[-0.5,1.2],'FontSize',ABCLabelFontSize);
                else
                    ylabel('Hz');
                    xlabel('deg');
                    t = title('B','Units','normalized','Position',[-0.5,1.2],'FontSize',ABCLabelFontSize);
                end
            end
        end
    end
    saveas(gcf,'test002.fig');
end

if testResult
%% Results
%
%   single         avg.
% E
%   1         3    9 
%    4x2      4    7
% I     
%   2         5    10       
%    2x2      6    8


LineScheme = {':','-.','--','-',':'};
dtheta = 180/ntheta;
evenOdd = mod(ntheta,2);
relTheta = (-(ntheta+evenOdd)/2:(ntheta+evenOdd)/2)*dtheta;
hResult = figure;
evenOdd = mod(ntheta,2);

if ~test
%1 and 2
pos = {[1,2,7,8,13,14,19,20],[25,26,31,32]};
columns = 2;
nl(1) = length(pos{1});
nl(2) = length(pos{2});
rows = nl./columns;
for k = 1:2
    if k==1
        neuronlist = randi(p.nv1e,[nl(k),1]);
        level = 2;
        type0 = [true(p.nv1e,1);false(p.nv1i,1)];
        types = type0(PKsortedID(:,level));
        candidates = find(types);
        ncand = length(candidates);
        neuronlist(1) = PKsortedID(candidates(1),level);
        neuronlist(3) = PKsortedID(candidates(floor(ncand/2)),level);
        neuronlist(5) = PKsortedID(candidates(ncand),level);
        level = 1;
        efired = nP(level).pkrate(dCVsortedID) > nP(level).br(dCVsortedID) & nP(level).ei(dCVsortedID) > 0.5 &...
                 nP(level).pkrate(dCVsortedID) > thres;
        type0 = [true(p.nv1e,1); false(p.nv1i,1)];
        types = type0(dCVsortedID);
        candidates = find(efired(dCVsortedID) & types); 
        ncand = length(candidates); 
        neuronlist(2) = dCVsortedID(candidates(2));
        neuronlist(4) = dCVsortedID(candidates(floor(ncand/2)));
        neuronlist(6) = dCVsortedID(candidates(ncand-1));

        level = 1;
        type0 = [true(p.nv1e,1);false(p.nv1i,1)];
        types = type0(PKsortedID(:,level));
        candidates = find(types);
        ncand = length(candidates);
        neuronlist(7) = PKsortedID(candidates(ncand),level);
        neuronlist(8) = PKsortedID(candidates(1),level);
        disp(neuronlist);
    else
        level = 2;
        ifired = nP(level).pkrate(dCVsortedID) > nP(level).br(dCVsortedID) & nP(level).ei(dCVsortedID) < 0.5 &...
                 nP(level).pkrate(dCVsortedID) > thres;
        neuronlist = p.nv1e + randi(p.nv1i,[nl(k),1]);
        type0 = [false(p.nv1e,1); true(p.nv1i,1)];
        types = type0(dCVsortedID);
        candidates = find(ifired(dCVsortedID) & types); 
        ncand = length(candidates); 
        neuronlist(1) = dCVsortedID(candidates(floor(ncand/1.2)));
        neuronlist(2) = dCVsortedID(candidates(floor(ncand/2)));
        neuronlist(3) = dCVsortedID(candidates(ncand));
        disp(neuronlist);
    end
    for j = 1:nl(k)
        i = neuronlist(j); 
        ax = subplot(sum(rows),columns*3,pos{k}(j));
        hold on
        inputAngle = nP(contrastLevel).priA(i);
        ipA = round(inputAngle*180/pi/dtheta)+1;
        if ipA < ntheta/2+1
            half = [(ipA+ntheta/2):ntheta,1:(ipA+ntheta/2)];
        else 
            half = [(ipA-ntheta/2):ntheta,1:(ipA-ntheta/2)];
        end
        if inputAngle >= pi/2
            pA = inputAngle - pi/2;
        else
            pA = inputAngle + pi/2;
        end
        grating_ipA = round(pA*180/pi/dtheta)+1;
        thetas = ((grating_ipA-ntheta/2-1):(grating_ipA+ntheta/2-1))*dtheta;
        for ci=1:contrastLevel
    
            rho = [tC(ci).frate(i,:),tC(ci).frate(i,1)];
            plot(thetas,rho(half),'Color','k','LineStyle',LineScheme{ci});
        end
        xlim([thetas(1),thetas(end)]);
        ylim([0,1.1*max(rho)]);
        set(ax,'XTick',[thetas(1),thetas(1)+90,thetas(end)]);
     
        if j == 1 
            if k == 1
                ylabel('Hz');
                xlabel('deg');
                t = title('A','Units','normalized','Position',[-0.5,1.2],'FontSize',ABCLabelFontSize);
            else
                ylabel('Hz');
                xlabel('deg');
                t = title('B','Units','normalized','Position',[-0.5,1.2],'FontSize',ABCLabelFontSize);
            end
        end
    end
end

%3 avg E
rasterTC = zeros(ntheta+1,p.nv1,contrastLevel)-1;
EmeanTC = zeros(ntheta+1,contrastLevel)-1;
EquaLTC = zeros(ntheta+1,contrastLevel)-1;
EquaUTC = zeros(ntheta+1,contrastLevel)-1;
EmediTC = zeros(ntheta+1,contrastLevel)-1;
ImeanTC = zeros(ntheta+1,contrastLevel)-1;
IquaLTC = zeros(ntheta+1,contrastLevel)-1;
IquaUTC = zeros(ntheta+1,contrastLevel)-1;
ImediTC = zeros(ntheta+1,contrastLevel)-1;
TCgLGN = rasterTC;
TCgLGN_F1 = rasterTC;
TCgE = rasterTC;
TCgI = rasterTC;
EmeanTCgLGN = EmeanTC;
ImeanTCgLGN = ImeanTC;
EmeanTCgLGN_F1 = EmeanTC;
ImeanTCgLGN_F1 = ImeanTC;
EmeanTCgE = EmeanTC;
ImeanTCgE = ImeanTC;
EmeanTCgI = EmeanTC;
ImeanTCgI = ImeanTC;
for i = 1:contrastLevel
    epick = nP(contrastLevel).ei>0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    ipick = nP(contrastLevel).ei<0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    for k=1:p.nv1    
        ipA = round(nP(i).priA(k)*180/pi/dtheta)+1;
        %ipA = round(nP(i).prA(k)*180/pi/dtheta)+1;
        if ipA > ntheta
            ipA = ipA-ntheta;
        end
        if ipA < ntheta/2+1
            half = [(ipA+ntheta/2):ntheta,1:ipA+ntheta/2];
        else 
            half = [ipA-ntheta/2:ntheta,1:(ipA-ntheta/2)];
        end
        rho = [tC(i).frate(k,:),tC(i).frate(k,1)];
        gLGN_F1 = [tC(i).gLGN(k,:).*tC(i).gLGNsc(k,:), tC(i).gLGN(k,1)*tC(i).gLGNsc(k,1)];
        gLGN = [tC(i).gLGN(k,:), tC(i).gLGN(k,1)];
        gE = [tC(i).gE(k,:),tC(i).gE(k,1)];
        gI = [tC(i).gI(k,:),tC(i).gI(k,1)];
        rasterTC(:,k,i) = rho(half);
        TCgLGN(:,k,i) = gLGN(half);
        TCgLGN_F1(:,k,i) = gLGN_F1(half);
        TCgE(:,k,i) = gE(half);
        TCgI(:,k,i) = gI(half);
    end
    EmeanTC(:,i) = mean(rasterTC(:,epick,i),2);
    EmediTC(:,i) = median(rasterTC(:,epick,i),2);
    EquaLTC(:,i) = EmediTC(:,i)-quantile(rasterTC(:,epick,i),0.25,2);
    EquaUTC(:,i) = quantile(rasterTC(:,epick,i),0.75,2)-EmediTC(:,i);
    EmeanTCgLGN(:,i) = mean(TCgLGN(:,epick,i),2);
    EmeanTCgLGN_F1(:,i) = mean(TCgLGN_F1(:,epick,i),2);
    EmeanTCgE(:,i) = mean(TCgE(:,epick,i),2);
    EmeanTCgI(:,i) = mean(TCgI(:,epick,i),2);

    ImeanTC(:,i) = mean(rasterTC(:,ipick,i),2);
    ImediTC(:,i) = median(rasterTC(:,ipick,i),2);
    IquaLTC(:,i) = ImediTC(:,i)-quantile(rasterTC(:,ipick,i),0.25,2);
    IquaUTC(:,i) = quantile(rasterTC(:,ipick,i),0.75,2) - ImediTC(:,i);
    ImeanTCgLGN(:,i) = mean(TCgLGN(:,ipick,i),2);
    ImeanTCgLGN_F1(:,i) = mean(TCgLGN_F1(:,ipick,i),2);
    ImeanTCgE(:,i) = mean(TCgE(:,ipick,i),2);
    ImeanTCgI(:,i) = mean(TCgI(:,ipick,i),2);
end

ax1 = subplot(4,3,2);
    ip = get(ax1,'Position');
    set(ax1,'Position',[ip(1)+0.03,ip(2)+0.02,ip(3)-0.09,ip(4)]);
    ip = get(ax1,'Position');
    ax2 = axes('Position',ip,'YAxisLocation','right','YColor','b');
    set(ax1,'Color','None');
    for i=1:contrastLevel
        axes(ax1);
        hold(ax1,'on')
        plot(relTheta, EmeanTC(:,i),LineScheme{i},'Color','k');
        axes(ax2);
        hold(ax2,'on');
        if i == contrastLevel
            %h(1) = plot(relTheta, EmeanTCgLGN(:,i),LineScheme{i},'Color','g');
            h(1) = plot(relTheta, EmeanTCgLGN_F1(:,i),LineScheme{i},'Color','m');
            h(2) = plot(relTheta, EmeanTCgE(:,i),LineScheme{i},'Color','r');
        else
            %plot(relTheta, EmeanTCgLGN(:,i),LineScheme{i},'Color','g');
            plot(relTheta, EmeanTCgLGN_F1(:,i),LineScheme{i},'Color','m');
            plot(relTheta, EmeanTCgE(:,i),LineScheme{i},'Color','r');
        end
        %plot(relTheta, EmeanTCgI(:,i),LineScheme{i},'Color','b');
    end
    %xlabel(ax1,'relative angle to pref(^{o})');
    xlabel(ax1,'deg');

    %legend(ax2,h,{'g_{LGN}','gLGN-F_1','g_E'},'Location','NorthEast');
    ylabel(ax1,'Hz'); ylabel(ax2,'Conductance');
    %yy = ylim(ax1);
    %ylim(ax1,[0,1.0*yy(2)]);
    %yy = ylim(ax2);
    %ylim(ax2,[0,1.0*yy(2)]);
    axes(ax1);
    xticks = -90:45:90;
    xlim(ax1,[relTheta(1),relTheta(end)]);
    xlim(ax2,[relTheta(1),relTheta(end)]);
    set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
    set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
    t = title('C','Units','normalized','Position',[-0.3,1.05],'FontSize',ABCLabelFontSize);

%4 norm E
ax = subplot(4,3,5);
    ip = get(ax,'Position');
    set(ax,'Position',[ip(1)+0.03,ip(2)+0.02,ip(3)-0.09,ip(4)]);
    hold on
    maxEmeanTC = ones(ntheta+1,1) * max(EmeanTC);
    maxEmeanTCgLGN = ones(ntheta+1,1) * max(EmeanTCgLGN);
    maxEmeanTCgLGN_F1 = ones(ntheta+1,1) * max(EmeanTCgLGN_F1);
    maxEmeanTCgE = ones(ntheta+1,1) * max(EmeanTCgE);
    maxEmeanTCgI = ones(ntheta+1,1) * max(EmeanTCgI);
    for i = 1:contrastLevel
        %plot(relTheta, EmeanTCgLGN(:,i)./maxEmeanTCgLGN(:,i),LineScheme{i},'Color','g');
        plot(relTheta, EmeanTCgLGN_F1(:,i)./maxEmeanTCgLGN(:,i),LineScheme{i},'Color','m');
        plot(relTheta, EmeanTCgE(:,i)./maxEmeanTCgE(:,i),LineScheme{i},'Color','r');
        plot(relTheta, EmeanTC(:,i)./maxEmeanTC(:,i),LineScheme{i},'Color','k');
    end
    ylabel('normalized');
    xlabel('deg');
    yy = ylim();
    ylim([0,yy(2)]);
    xticks = -90:45:90;
    xlim([relTheta(1),relTheta(end)]);
    set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'));
    t = title('D','Units','normalized','Position',[-0.3,1.05],'FontSize',ABCLabelFontSize);

%5 avg I
ax1 = subplot(4,3,8);
    ip = get(ax1,'Position');
    set(ax1,'Position',[ip(1)+0.03,ip(2)+0.02,ip(3)-0.09,ip(4)]);
    ip = get(ax1,'Position');
    ax2 = axes('Position',ip,'YAxisLocation','right','YColor','b');
    axes(ax1);
    for i=1:contrastLevel
        hold(ax1,'on')
        plot(relTheta, ImeanTC(:,i),LineScheme{i},'Color','k');
    end
    axes(ax2);
    hold(ax2,'on');
    for i=1:contrastLevel
        %plot(relTheta, ImeanTCgLGN(:,i),LineScheme{i},'Color','g');
        plot(relTheta, ImeanTCgLGN_F1(:,i),LineScheme{i},'Color','m');
        plot(relTheta, ImeanTCgE(:,i),LineScheme{i},'Color','r');
        %plot(relTheta, EmeanTCgI(:,i),LineScheme{i},'Color','b');
    end
    set(ax2,'YColor','b')
    ylabel(ax1,'Hz'); ylabel(ax2,'Conductance');
    xlabel(ax1,'deg');
    yy = ylim(ax1);
    ylim(ax1,[0,yy(2)]);
    yy = ylim(ax2);
    ylim(ax2,[0,yy(2)]);
    axes(ax1);
    xticks = -90:45:90;
    xlim(ax1,[relTheta(1),relTheta(end)]);
    xlim(ax2,[relTheta(1),relTheta(end)]);
    set(ax1,'XTick',xticks,'XTickLabel',num2str(xticks'));
    set(ax2,'XTick',xticks,'XTickLabel',num2str(xticks'));
    %set(ax2,'XTick',[],'XTickLabel','','XTickLabelMode','manual');
    set(ax1,'Color','None');
    t = title('E','Units','normalized','Position',[-0.3,1.05],'FontSize',ABCLabelFontSize);

%6 norm I
ax = subplot(4,3,11);
    ip = get(ax,'Position');
    set(ax,'Position',[ip(1)+0.03,ip(2)+0.02,ip(3)-0.09,ip(4)]);
    hold on
    maxImeanTC = ones(ntheta+1,1) * max(ImeanTC);
    maxImeanTCgLGN = ones(ntheta+1,1) * max(ImeanTCgLGN);
    maxImeanTCgLGN_F1 = ones(ntheta+1,1) * max(ImeanTCgLGN_F1);
    maxImeanTCgE = ones(ntheta+1,1) * max(ImeanTCgE);
    maxImeanTCgI = ones(ntheta+1,1) * max(ImeanTCgI);
    for i = 1:contrastLevel
        %plot(relTheta, ImeanTCgLGN(:,i)./maxImeanTCgLGN(:,i),LineScheme{i},'Color','g');
        plot(relTheta, ImeanTCgLGN_F1(:,i)./maxImeanTCgLGN(:,i),LineScheme{i},'Color','m');
        plot(relTheta, ImeanTCgE(:,i)./maxImeanTCgE(:,i),LineScheme{i},'Color','r');
        plot(relTheta, ImeanTC(:,i)./maxImeanTC(:,i),LineScheme{i},'Color','k');
    end
    ylabel('norm. max');
    xlabel('deg');
    yy = ylim();
    ylim([0,yy(2)]);
    title('Inh');
    xticks = -90:45:90;
    xlim([relTheta(1),relTheta(end)]);
    set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'));
    t = title('F','Units','normalized','Position',[-0.3,1.05],'FontSize',ABCLabelFontSize);
end
%9 C2-C4 E
hEmax = 1;
hEmin = 0;
hImax = 0.4;
hImin = 0;
p1 =contrastLevel-2; p2 = contrastLevel;
while p1<=0
    p1 = p1 + 1;
end
nCV = 8;
ctrs = cell(2,1);
dTick = 0.2;
pick = nP(contrastLevel).ei>0.5 & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(p1).pkrate>thres & nP(p2).pkrate>thres;
pair1 = 1-nP(p1).cv(pick);
pair2 = 1-nP(p2).cv(pick);
%minCtrs = min(hEmin,min(min(pair1),min(pair2)));
%maxCtrs = max(hEmax,max(max(pair1),max(pair2)));
minCtrs = 0;
maxCtrs = 1;
%[tick,n0] = autoAxis(minCtrs,maxCtrs,round(1/dTick),[0,1]);
tick = 0:0.2:1.0;
n0 = length(tick);
tickLabel = num2str(tick');
lctrsx = (n0-1)*nCV+1;
lctrsy = lctrsx;
ctrs{1} = linspace(tick(1),tick(end),lctrsx);
ctrs{2} = ctrs{1};
tickPosY = 0.5:nCV:(lctrsy-1+0.5);
tickPosX = 0.5:nCV:(lctrsx-1+0.5);
hExc = subplot(4,3,3);
    ep = get(hExc,'Position');
    set(hExc,'Position',[ep(1)-0.01,ep(2)-0.05,ep(3),ep(3)]);
	CVpair = [pair1, pair2];
    denCVpair = hist3(CVpair,ctrs);
    maxDen = max(max(denCVpair));
    denCVpair = denCVpair/maxDen;
    denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
    img = imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
    hold on
    plot(lctrsx:-1:0,0:lctrsy,'-.k','LineWidth',2);
    %axis(hExc,'ij');
    %set(hExc,'Box','on');
    %axis(hExc,'image');
    %axis(hExc,'normal');
    xlabel('gOSI(25%)')
    ylabel('gOSI(100%)')
    colormap(hExc,redOnly);
    
    set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
    t = title('G','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);
%end
dTick = 0.1;
nCV = 10;
%10 C2-C4 I
pick = nP(contrastLevel).ei<0.5 & nP(p1).pkrate> nP(1).br & nP(p2).pkrate> nP(1).br & nP(p1).pkrate>thres & nP(p2).pkrate>thres;
pair1 = 1-nP(p1).cv(pick);
pair2 = 1-nP(p2).cv(pick);
minCtrs = max(hImin,min(min(pair1),min(pair2)));
maxCtrs = min(hImax,max(max(pair1),max(pair2)));
%[tick,n0] = autoAxis(minCtrs,maxCtrs,round(hImax/dTick),[0,1]);
tick = 0:0.1:0.4;
n0 = length(tick);
%disp(tick);
tickLabel = num2str(tick');
lctrsx = (n0-1)*nCV+1;
lctrsy = lctrsx;
ctrs{1} = linspace(tick(1),tick(end),lctrsx);
ctrs{2} = ctrs{1};
%tickPosY =(tick-tick(1))./(tick(end)-tick(1))*lctrsy;
%tickPosX =(tick-tick(1))./(tick(end)-tick(1))*lctrsx;
tickPosY = 0.5:nCV:(lctrsx-1+0.5);
tickPosX = 0.5:nCV:(lctrsy-1+0.5);
hInh = subplot(4,3,9);
    ip = get(hInh,'Position');
    set(hInh,'Position',[ip(1)-0.01,ip(2)-0.05,ip(3),ip(3)]);
    
    CVpair = [pair1, pair2];
    denCVpair = hist3(CVpair,ctrs);
    maxDen = max(max(denCVpair));
    denCVpair = denCVpair/maxDen;
    denCVpair = denCVpair(1:(lctrsx-1),1:(lctrsy-1));
    img = imagesc([1,lctrsx-1],[lctrsy-1,1],denCVpair');
    %cd = rand(lctrsx-1,lctrsy-1);
    %cd = triu(cd);
    %img = imagesc([1,lctrsx-1],[lctrsy-1,1],cd);
    hold on
    plot(lctrsx:-1:0,0:lctrsy,'-.k','LineWidth',2);
    %axis(hInh,'ij');
    %axis(hInh,'image');
    %axis(hInh,'normal');
    xlabel('gOSI(25%)')
    ylabel('gOSI(100%)')
    colormap(hInh,blueOnly);
    
    set(gca,'YTickLabel',flipud(tickLabel),'YTick',tickPosY,'XTickLabel',tickLabel,'XTick',tickPosX);
    t = title('I','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);
    %save([opfdr,'/tempRedraw.mat'],'denCVpair','lctrsx','lctrsy');
if ~test
% Exc F1/F0
ranges = 0:0.2:2; 
xSimp = 1.1:0.2:1.9;
xComp = 0.1:0.2:0.9;
ax = subplot(4,3,6);
    ip = get(ax,'Position');
    set(ax,'Position',[ip(1)-0.01,ip(2),ip(3),ip(4)*0.7]);
    hold on
    counts = histcounts(nP(contrastLevel).sc(1:p.nv1e),ranges);
    counts = counts./sum(counts)*100;
    bar(xSimp,counts(6:10),'w');
    bar(xComp,counts(1:5),'k');
    ylabel('%');
    xlabel('F1/F0');
    t = title('H','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);
% Inh F1/F0
ax = subplot(4,3,12);
    ip = get(ax,'Position');
    set(ax,'Position',[ip(1)-0.01,ip(2),ip(3),ip(4)*0.7]);
    hold on
    counts = histcounts(nP(contrastLevel).sc(p.nv1e+(1:p.nv1i)),ranges);
    counts = counts./sum(counts)*100;
    bar(xSimp,counts(6:10),'w');
    bar(xComp,counts(1:5),'k');
    ylabel('%');
    xlabel('F1/F0');
    t = title('J','Units','normalized','Position',[-0.2,1.05],'FontSize',ABCLabelFontSize);
    %set(ax,'Visible','off')
end
if ~isempty(format)
    %get(gcf,'renderer')
    set(gcf,'renderer','painters');
    set(gcf,'PaperUnits','points');
    set(gcf,'PaperSize',[ResultW,ResultH]);
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperPosition',[0,0,ResultW,ResultH]);
    set(gcf,'Units','points');
    set(gcf,'OuterPosition',[ResultMW0,ResultMH0,ResultW-ResultMW1-ResultMW0,ResultH-ResultMH1-ResultMH0]);
    print(hResult,[opfdr,'/Results-',theme,'.',format],printDriver,dpi);
end
saveas(hResult,[opfdr,'/Results-',theme,'.fig']);
end
%%
if test000
    LineScheme = {':','-.','--','-',':'};
    dtheta = 180/ntheta;
    evenOdd = mod(ntheta,2);
    relTheta = (-(ntheta+evenOdd)/2:(ntheta+evenOdd)/2)*dtheta;
    load(data0);
    EmeanTCgI = zeros(ntheta+1,contrastLevel)-1;
TCgI = zeros(ntheta+1,p.nv1,contrastLevel)-1;
for i = 1:contrastLevel
    epick = nP(contrastLevel).ei>0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    for k=1:p.nv1    
        ipA = round(nP(i).priA(k)*180/pi/dtheta)+1;
        %ipA = round(nP(i).prA(k)*180/pi/dtheta)+1;
        if ipA > ntheta
            ipA = ipA-ntheta;
        end
        if ipA < ntheta/2+1
            half = [(ipA+ntheta/2):ntheta,1:ipA+ntheta/2];
        else 
            half = [ipA-ntheta/2:ntheta,1:(ipA-ntheta/2)];
        end
        gI = [tC(i).gI(k,:),tC(i).gI(k,1)];
        TCgI(:,k,i) = gI(half);
    end
    EmeanTCgI(:,i) = mean(TCgI(:,epick,i),2);
end
    
load(dataf);
EmeanTCgIf = zeros(ntheta+1,contrastLevel)-1;
TCgIf = zeros(ntheta+1,p.nv1,contrastLevel)-1;
for i = 1:contrastLevel
    epick = nP(contrastLevel).ei>0.5 & nP(i).pkrate> nP(1).br & nP(contrastLevel).pkrate>thres;
    for k=1:p.nv1    
        ipA = round(nP(i).priA(k)*180/pi/dtheta)+1;
        %ipA = round(nP(i).prA(k)*180/pi/dtheta)+1;
        if ipA > ntheta
            ipA = ipA-ntheta;
        end
        if ipA < ntheta/2+1
            half = [(ipA+ntheta/2):ntheta,1:ipA+ntheta/2];
        else 
            half = [ipA-ntheta/2:ntheta,1:(ipA-ntheta/2)];
        end
        gI = [tC(i).gI(k,:),tC(i).gI(k,1)];
        TCgIf(:,k,i) = gI(half);
    end
    EmeanTCgIf(:,i) = mean(TCgIf(:,epick,i),2);
end

figure;
axgI = subplot(1,2,1);
%     ip = get(ax,'Position');
    %plot(1,1,'*');
    i=contrastLevel;
    %objectMax = max(EmeanTCgI(:,i));
    hold on
    for j=1:contrastLevel
    %for j=1:contrastLevel-1
        %plot(relTheta, EmeanTCgI(:,j)*objectMax/max(EmeanTCgI(:,j)),LineScheme{j},'Color','b');
        plot(relTheta, EmeanTCgI(:,j)./max(EmeanTCgI(:,j)),LineScheme{j},'Color','b');
        plot(relTheta, EmeanTCgIf(:,j)./max(EmeanTCgIf(:,j)),LineScheme{j},'Color','k');
    end
    %plot(relTheta, EmeanTCgI(:,i),LineScheme{i},'Color','b');
    xticks = -90:45:90;
    xlim([relTheta(1),relTheta(end)]);
    xlabel('deg','FontSize',8);
    set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'),'FontSize',FontSize);
%     set(ax,'Position',[ip(1)+0.03,ip(2)+0.01,ip(3)*0.9,ip(4)*0.9]);
    %set(ax,'Position',[ip(1),ip(2),ip(3),ip(4)]);
    ylabel('normalized g_{I}');
    %legend(conLabel);
    %title(['Normalized to C',conLabel{i}]);
    %t = title(ax,'F','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);
    subplot(1,2,2);
    hold on
    for j=1:contrastLevel
    %for j=1:contrastLevel-1
        %plot(relTheta, EmeanTCgI(:,j)*objectMax/max(EmeanTCgI(:,j)),LineScheme{j},'Color','b');
        yyaxis left
        plot(relTheta, EmeanTCgI(:,j)./max(EmeanTCgI(:,j)),LineScheme{j},'Color','b');
        yyaxis right
        plot(relTheta, EmeanTCgIf(:,j),LineScheme{j},'Color','k');
    end
    xticks = -90:45:90;
    xlim([relTheta(1),relTheta(end)]);
    xlabel('deg','FontSize',8);
    set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'),'FontSize',FontSize);
end
%%
if testAnalysis


%Analysis
%
%       2   5
%           
%   1   3   6
%            
%       4   7
hAnalysis = figure;

if rr == 4
    fdr0 = {'isample40-tf1l','isample40-t101l','isample40-t81l','isample40-t61l'};
    fdr1 = {'isample40-tf2el','isample40-t102el','isample40-t82el','isample40-t62el'};
    fdr2 = {'isample40-nf1l','isample40-n101l','isample40-n81l','isample40-n61l'};
end
if rr == 5
    fdr0 = {'sample-5tf5021','sample-5t5021','sample-5t8021','sample-5t6021'};
    fdr1 = {'sample-5tf1-ie','sample-5t101-ie','sample-5t81-ie','sample-5t61-ie'};
    fdr2 = {'sample-5tfe1','sample-5te1','sample-5t8e1','sample-5t6e1'};
end
nfdr = length(fdr0);

pLabel = '\sigma of I{\rightarrow}E';
pTickLabel = {'\sigma_{I{\rightarrow}E} = \infty','\sigma_{I{\rightarrow}E} = 1.0','\sigma_{I{\rightarrow}E} = 0.8','\sigma_{I{\rightarrow}E} = 0.6'};
mks = 26;
%exp
tt0 = true 
if tt0 
hExp = subplot(3,6,1);
    ip = get(hExp,'Position');
    set(hExp,'Position',[ip(1)-0.0,ip(2)+0.068,ip(3)*0.5,ip(4)*0.54]);
    img = imread('sigma.png');
    image(img);
    axis(hExp,'off');
    t = title('A','Units','normalized','Position',[-1.0,1.00],'FontSize',ABCLabelFontSize);
%1
hMech = subplot(3,3,[4,7]);
    ip = get(hMech,'Position');
    set(hMech,'Position',[0.02,ip(2)-0.05,ip(3)*1.7,ip(4)*1.3]);
    img = imread('mechanism.png');
    image(img);
    axis(hMech,'off');
    t = title('B','Units','normalized','Position',[0.10,0.90],'FontSize',ABCLabelFontSize);
end

%2
tt1 = true 
if tt1
ax = subplot(3,3,2);
    ip = get(ax,'Position');
    delete(ax);
    [hCV0, hCond0, hFR0, hFRr0, hgCV0] = parameterTrace(fdr0,pLabel,pTickLabel,lgnfile,'','',mks);
    hCV = findobj(hCV0,'Type','axes');
    hCVE = hCV(2);
    leg = findobj(hCV0,'Type','legend');
    legE = leg(2);
    handles = copyobj([hCVE,legE],hAnalysis);
    ax = handles(1);
    lg = handles(2);
    set(ax,'Position',[ip(1)+0.04,ip(2),ip(3)*0.85,ip(4)]);
    set(ax,'FontSize',8);
    set(lg,'FontSize',8);
    lg.Box = 'off';
    ylim(ax,[0,1]);
    xlim(ax,[0,1]);
    lp = get(lg,'Position');
    set(lg,'Position',[lp(1)-0.34,lp(2)+0.07,lp(3)*0.85,lp(4)]);
    t = title(ax,'C','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);
    txt = text('Position',[0.1,1.0,0],'String','Standard','Units','normalized','FontSize',8);
    txt.Parent = ax;
    set(ax,'YTick',[0,0.5,1]),
%3
ax = subplot(3,3,5);
    ip = get(ax,'Position');
    delete(ax);
    [hCV1, hCond1, ~,~] = parameterTrace(fdr1,pLabel,pTickLabel,lgnfile,'','',mks);
    hCV = findobj(hCV1,'Type','axes');
    hCVE = hCV(2);
    handles = copyobj(hCVE,hAnalysis);
    ax = handles(1);
    set(ax,'Position',[ip(1)+0.04,ip(2),ip(3)*0.85,ip(4)]);
    set(ax,'FontSize',8);
    ylim(ax,[0,1]);
    xlim(ax,[0,1]);
    t = title(ax,'D','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);
    txt = text('Position',[0.1,0.9,0],'String','70% inhibtion','Units','normalized','FontSize',8);
    txt.Parent = ax;
    set(ax,'YTick',[0,0.5,1]),

%4
ax = subplot(3,3,8);
    ip = get(ax,'Position');
    delete(ax);
    [hCV2, hCond2, ~,~] = parameterTrace(fdr2,pLabel,pTickLabel,lgnfile,'','',mks);
    hCV = findobj(hCV2,'Type','axes');
    hCVE = hCV(2);
    handles = copyobj(hCVE,hAnalysis);
    ax = handles(1);
    set(ax,'Position',[ip(1)+0.04,ip(2),ip(3)*0.85,ip(4)]);
    set(ax,'FontSize',8);
    ylim(ax,[0,1]);
    xlim(ax,[0,1]);
    t = title(ax,'E','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);
    txt = text('Position',[0.1,0.9,0],'String','Uniform E{\leftarrow}E','Units','normalized','FontSize',8);
    txt.Parent = ax;
    set(ax,'YTick',[0,0.5,1]),
else
    [hCV0, hCond0, hFR0, hFRr0] = parameterTrace(fdr0,pLabel,pTickLabel,lgnfile,'','',mks);
    [hCV1, hCond1, ~,~] = parameterTrace(fdr1,pLabel,pTickLabel,lgnfile,'','',mks);
    [hCV2, hCond2, ~,~] = parameterTrace(fdr2,pLabel,pTickLabel,lgnfile,'','',mks);
end 

tt = true
if tt
%6
ax = subplot(4,3,6);
    ip = get(ax,'Position');
    delete(ax);
    hFRrs = findobj(hFRr0,'Type','axe');
    hFRr = hFRrs(1);
    handles = copyobj(hFRr,hAnalysis);
    ax = handles(1);
    set(ax,'Position',[ip(1)+0.03,ip(2),ip(3)*0.9,ip(4)]);
    %set(ax,'Position',[ip(1),ip(2),ip(3),ip(4)]);
    set(ax,'FontSize',8);
    t = title(ax,'G','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);

%7
ax = subplot(4,3,9);
    ip = get(ax,'Position');
    delete(ax);
    hCond = findobj(hCond0,'Type','axes');
    %lCond = findobj(hCond0,'Type','legend');
    hCondE = hCond(2);
    %lCondE = lCond(2);
    %handles = copyobj([hCondE,lCondE],hAnalysis);
    handles = copyobj(hCondE,hAnalysis);
    ax = handles(1);
    lines = findobj(ax,'Type','line');

    set(ax,'Position',[ip(1)+0.03,ip(2),ip(3)*0.9,ip(4)]);
    ylim(ax,[0.2,1.21]);
    %set(ax,'Position',[ip(1),ip(2),ip(3),ip(4)]);
    %lg = handles(2);
    %lg = legend(ax,lines(2:-1:1),{'C 25%','C 100%'},'Location','NorthEast','FontSize',8,'Units','normalized');
    %lp = get(lg,'Position');
    %set(lg,'Position',[lp(1)-0.015,lp(2)+0.01,lp(3),lp(4)]);
    set(ax,'FontSize',8);
    t = title(ax,'H','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);

%8
ax = subplot(4,3,12);
    ip = get(ax,'Position');
    delete(ax);
    hFRs = findobj(hFR0,'Type','axe');
    hFR = hFRs(1);
    handles = copyobj(hFR,hAnalysis);
    ax = handles(1);
    set(ax,'Position',[ip(1)+0.03,ip(2),ip(3)*0.9,ip(4)]);
    %set(ax,'Position',[ip(1),ip(2),ip(3),ip(4)]);
    set(ax,'FontSize',8);
    t = title(ax,'I','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);
end

%5
%ax = subplot(4,3,3);
%    ip = get(ax,'Position');
%    delete(ax);
%    hgCV = findobj(hgCV0,'Type','axes');
%    hgCVI = hgCV(1);
%    handles = copyobj(hgCVI,hAnalysis);
%    ax = handles(1);
%    set(ax,'FontSize',8);
%    set(ax,'Position',[ip(1)+0.03,ip(2),ip(3)*0.9,ip(4)*0.9]);
%    t = title(ax,'F','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);
%    ylim(ax,[0,0.1]);
%    xlim(ax,[0,0.1]);
%    txt = text('Position',[0.1,0.9,0],'String','inhibitory conductance','Units','normalized','FontSize',8);
%    txt.Parent = ax;
%    set(ax,'YTick',[0,0.5,1]),
%    t = title(ax,'F','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);

clear hFR0 hCV0 hCond0 hFRr0 hCV1 hCond1 hCV2 hCond2
ax = subplot(4,3,3);
    ip = get(ax,'Position');
    delete(ax);
    ax = copyobj(axgI,hAnalysis);
    set(ax,'FontSize',8);
    set(ax,'Position',[ip(1)+0.03,ip(2),ip(3)*0.9,ip(4)*0.9]);
    t = title(ax,'F','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);

%ax = subplot(4,3,3);
%    ip = get(ax,'Position');
%    %plot(1,1,'*');
%    i=contrastLevel;
%    %objectMax = max(EmeanTCgI(:,i));
%    hold on
%    for j=1:contrastLevel
%    %for j=1:contrastLevel-1
%        %plot(relTheta, EmeanTCgI(:,j)*objectMax/max(EmeanTCgI(:,j)),LineScheme{j},'Color','b');
%        plot(relTheta, EmeanTCgI(:,j)./max(EmeanTCgI(:,j)),LineScheme{j},'Color','b');
%    end
%    %plot(relTheta, EmeanTCgI(:,i),LineScheme{i},'Color','b');
%    xticks = -90:45:90;
%    xlim([relTheta(1),relTheta(end)]);
%    set(gca,'XTick',xticks,'XTickLabel',num2str(xticks'));
%    set(ax,'Position',[ip(1)+0.03,ip(2),ip(3)*0.9,ip(4)]);
%    %set(ax,'Position',[ip(1),ip(2),ip(3),ip(4)]);
%    ylabel('normalized');
%    %legend(conLabel);
%    %title(['Normalized to C',conLabel{i}]);
%    t = title(ax,'F','Units','normalized','Position',[-0.25,0.95],'FontSize',ABCLabelFontSize);
%
saveas(hAnalysis,[opfdr,'/Analysis-',theme,'.fig']);
if ~isempty(format)
    set(gcf,'PaperUnits','points');
    set(gcf,'PaperSize',[AnalysisW,AnalysisH]);
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperPosition',[0,0,AnalysisW,AnalysisH]);
    set(gcf,'Units','points');
    set(gcf,'OuterPosition',[AnalysisMW0,AnalysisMH0,AnalysisW-AnalysisMW1-AnalysisMW0,AnalysisH-AnalysisMH1-AnalysisMH0]);
    print(hAnalysis,[opfdr,'/Analysis-',theme,'.',format],printDriver,dpi);
end
end
