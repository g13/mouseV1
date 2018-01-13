opfdr = 'PprFig';
type = 't06';
theme = 'ndi305-40';
lgnfile = ['1xu-',theme,'-s911'];
coMat = ['coMa-',lgnfile];
conMat = ['conMat-',theme,'_',type,'.mat'];
format = 'png';
boxOnOff = 'on';

LGN.Aa = 15.88 * (180/pi)^2;
LGN.Ab = Aa*0.97;
rc = 5.61; %deg
rs = 16.98*rc/5.61; %deg
rc = rc * pi/180;
rs = rs * pi/180;
LGN.siga2 = rc^2/2;
LGN.sigb2 = rs^2/2;

load(lgnfile);
fid = fopen(conMat);
    m = fread(fid,[p.nv1,p.nv1],'int8');
fclose(fid);
pPosition = [0, 0, 1280, 720];
LineWidth = 2;
set(groot,'defaultLineLineWidth',LineWidth);
set(groot,'defaultErrorbarLineWidth',LineWidth);
FontSize = 14;
set(groot,'defaultAxesFontSize',FontSize);
set(groot,'defaultTextFontSize',FontSize);
LegendOffset = 1;
set(groot,'defaultLegendFontSize',FontSize-LegendOffset);
if ~isempty(format)
    if strcmp(format,'psc2')
        printDriver = ['-de',format];
        format = 'eps';
    else
        printDriver = ['-d',format];
    end
    dpi = '-r100';
end

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
hLGN = figure;

%   1.1 RF_E
LW = 5;
dxy = 24;
redge = 1;
i = randi(p.nv1e,1);
ax = subplot(2,2,1);
    hold on
    nSubregion = p.eSubregion(i);
    xlimit = [inf,-inf];
    ylimit = xlimit;
    for k = 1:nSubregion
        xtmp = [min(LGNpos(v1Map{i,k},1))-LGN.rs*redge, max(LGNpos(v1Map{i,k},1)) + LGN.rs*redge];
        if xtmp(1) < xlimit(1), xlimit(1) = xtmp(1); end
        if xtmp(2) > xlimit(2), xlimit(2) = xtmp(2); end
        ytmp = [min(LGNpos(v1Map{i,k},1))-LGN.rs*redge, max(LGNpos(v1Map{i,k},1)) + LGN.rs*redge];
        if ytmp(1) < ylimit(1), ylimit(1) = ytmp(1); end
        if ytmp(2) > ylimit(2), ylimit(2) = ytmp(2); end
    end
    x = xlimit(1):rs/dxy:xlimit(2);
    y = ylimit(1):rs/dxy:ylimit(2);
    for k = 1:nSubregion
        if sum(v1Map{k}) > 0
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(v1Map{i,k},:),x,y,nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','r','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','r','LineWidth',LW/2);
        else
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(-v1Map{i,k},:),x,y,-nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','b','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','b','LineWidth',LW/2);
        end
    end
    axis equal
    t = title('A','Units','normalized','Position', [0 1]);
    axis(ax,'off');

%   1.2 RF_I
i = p.nv1e + randi(p.nv1i,1);
ax = subplot(4,2,[2,4,6]);
    hold on
    nSubregion = p.iSubregion(i-p.nv1e);
    xlimit = [inf,-inf];
    ylimit = xlimit;
    for k = 1:nSubregion
        xtmp = [min(LGNpos(v1Map{i,k},1))-LGN.rs*redge, max(LGNpos(v1Map{i,k},1)) + LGN.rs*redge];
        if xtmp(1) < xlimit(1), xlimit(1) = xtmp(1); end
        if xtmp(2) > xlimit(2), xlimit(2) = xtmp(2); end
        ytmp = [min(LGNpos(v1Map{i,k},1))-LGN.rs*redge, max(LGNpos(v1Map{i,k},1)) + LGN.rs*redge];
        if ytmp(1) < ylimit(1), ylimit(1) = ytmp(1); end
        if ytmp(2) > ylimit(2), ylimit(2) = ytmp(2); end
    end
    x = xlimit(1):rs/dxy:xlimit(2);
    y = ylimit(1):rs/dxy:ylimit(2);
    for k = 1:nSubregion
        if sum(v1Map{k}) > 0
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(v1Map{i,k},:),x,y,nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','r','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','r','LineWidth',LW/2);
        else
            [X, Y, Z] = multiLGNspatialKernel_ext(LGNpos(-v1Map{i,k},:),x,y,-nSubLGN{i}(k),lgnStrength{i,k},LGN,200);
            contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','b','LineWidth',LW);
            newZmax = max(Z(:));
            contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','b','LineWidth',LW/2);
        end
    end
    axis equal
    axis(ax,'off');

%   1.3 Hist Norm. Distance
ax = subplot(2,2,3)
    md = max(p.enormDistance);
    gray = 0.5;
    dark = 0;
    edges = 0:0.1:ceil(md*10)/10;
    h = histogram(p.enormDistance, edges,'FaceColor',dark*ones(1,3),'EdgeColor',gray*ones(1,3));
    t = title('C','Units','normalized','Position', [0 1]);
    %ylabel('# Exc');
    xlabel('norm. distance');
    box(ax,'off');
    set(ax,'YColor','w','YTick',[]);

%   1.4 On-Off.png
subplot(4,6,7);
%   imread('On-Off.png');
    box(ax,'off');

fLGN = [opfdr,'/LGN-',theme];
saveas(hLGN,[fLGN,'.fig']);
set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
print(hLGN,[fLGN,'.png'],printDriver,dpi);

hSetup = figure;
%   1. LGN 
ax = subplot(2,3,1);
    imread([fLGN,'.png']);
    axis(ax,'off');
    
%   2. Salt-pepper (ellipses)
subplot(2,3,4);

%   3. Con. Prob. vs Distance 
subplot(3,3,2);

%   4. Con. Prob. vs Orientation Diff
subplot(3,3,5);

%   5. Con. Prob. vs RFcoeff 
subplot(3,3,8);

%   6. E-E logNormal 
subplot(3,3,3);

%   7. EPSP RFcoeff sample;
subplot(3,3,6);

%   8. Con. CDF, Exc. CDF vs RFcoeff
subplot(3,3,9);

saveas(hSetup,[opfdr,'/Setup-',theme,'.fig']);
if ~isempty(format)
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    print(hSetup,[opfdr,'/Setup-',theme,'.',format],printDriver,dpi);
end

