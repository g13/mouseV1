%testRun
function finalRF(lgnfile,eg,outputfdr,strengthFile,conMatFile)
if nargin < 4
    strengthFile = 'logNormalProfile.mat';
    if nargin < 5
        conMatFile = 'conMat.mat';
    end
end
% lgnfile = '1xu-n3t-s911.mat';
% conMatFile = 'conMat-1xu-n3t-s911-rf-gauss-coGauss-r50-s911n-logNormal.mat';
load(lgnfile);
load(strengthFile);
p2s = 20;
nLGNmean = mean(nLGN(1:p.nv1e,1));
matName = outputfdr;
fid = fopen(conMatFile);
m =fread(fid,[p.nv1,p.nv1],'int8');
nex = length(eg);
% nex = 10;
% rng(325);
% eg = [randi(p.nv1e,[nex,1]), randi(p.nv1i,[nex,1])];
% eg = [p.nv1e/2, p.nv1i/2+18];
format = 'fig';
format0 = 'png';
printDriver0 = ['-d',format0];
if ~isempty(format)
    if strcmp(format,'psc2')
        printDriver = ['-de',format];
        format = 'eps';
    else
        printDriver = ['-d',format];
    end
    dpi = '-r150';
end
pPosition = [18, 180, 1200, 900];
pex = ones(p.ev1y,1)*((1:p.ev1x)-0.5).*p.dE;
pex = reshape(pex,[p.nv1e,1]); % id goes along y first (column major)
pey = ((1:p.ev1y)-0.5)'*ones(1,p.ev1x).*p.dE;
pey = reshape(pey,[p.nv1e,1]);

pix = ones(p.iv1y,1) * (((1:p.iv1x)-0.5).*p.dI);
pix = reshape(pix,[p.nv1i,1]);
piy = (((1:p.iv1y)-0.5)'.*p.dI) * ones(1,p.iv1x);
piy = reshape(piy,[p.nv1i,1]);

elx = (p.ev1x) * p.dE;
ely = (p.ev1y) * p.dE;
ilx = (p.iv1x) * p.dI;
ily = (p.iv1y) * p.dI;




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
    i = eg(j);
    hExRF=figure;
    hold on
    subregion = length(nSubLGN{i});
    subpick = 1:subregion;
    if i <= p.nv1e
        theta = etheta(i);
		sigma = p.esigma(i,subpick);
    	sigmb = sigma.*p.eAspectRatio(i,subpick);
        EIlabel = 'E';
    else
        theta = itheta(i-p.nv1e);
		sigma = p.isigma(i-p.nv1e,subpick);
    	sigmb = sigma.*p.iAspectRatio(i-p.nv1e,subpick);
        EIlabel = 'I';
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
    edge = 4;
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
    contour(X,Y,Z_Eon,[0,0],'LineStyle',':','LineColor','r','LineWidth',2);
    contour(X,Y,Z_Eoff,[0,0],'LineStyle',':','LineColor','b','LineWidth',2);
    Z_on = Z_Eon + Z_Ion + Z_ff_on;
    Z_off = Z_Eoff + Z_Ioff + Z_ff_off;
    zcmax = max(max(max(Z_on)),max(max(Z_off)));
    
    mesh(X,Y,Z_on-Z_off,'FaceColor','none','EdgeColor','k');
    mesh(X,Y,Z_ff_on-Z_ff_off,'FaceColor','none','EdgeColor','g');
    
    hIon = mesh(X,Y,Z_Ion,'FaceColor','none','EdgeColor','r');
    set(hIon,'Visible','off');
    hIoff = mesh(X,Y,Z_Ioff,'FaceColor','none','EdgeColor','b');
    set(hIoff,'Visible','off');
    
    hEon = mesh(X,Y,Z_Eon,'FaceColor','none','EdgeColor','r');
    set(hEon,'Visible','off');
    hEoff = mesh(X,Y,Z_Eoff,'FaceColor','none','EdgeColor','b');
    set(hEoff,'Visible','off');
    
    hffon = mesh(X,Y,Z_ff_on,'FaceColor','none','EdgeColor','r');
    set(hffon,'Visible','off');
    hffoff = mesh(X,Y,Z_ff_off,'FaceColor','none','EdgeColor','b');
    set(hffoff,'Visible','off');
    
    [maxtmp, indi] = max(abs(Z_on-Z_off));
    [deltamax, indj] = max(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'*k');
    title({['Ex.',EIlabel,num2str(j)],['LGN RF max = ',num2str(zmax,'%3.1f'),', cortical RF max = ', num2str(zcmax,'%3.1f'),', max\delta=',num2str(deltamax,'%3.1f')]});
    
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
    
    [maxtmp, indi] = min(Z_Ion);
    [deltamax, indj] = min(maxtmp);
    plot3(x(indj),y(indi(indj)),deltamax,'vr');
    
    [maxtmp, indi] = min(Z_Ioff);
    [deltamax, indj] = min(maxtmp);
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
        figName = [outputfdr,'/',num2str(i),'-RF-',matName,'.',format];
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if ~strcmp(format,'fig')
            print(hExRF,figName,printDriver,dpi);
        else
            saveas(hExRF,figName);
        end
    end
end
mei = m(p.nv1e+(1:p.nv1i),1:p.nv1e);
for j = 1:nex
    hEx = figure;
    i = eg(j);
    if i <= p.nv1e
        subplot(1,2,1);
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
        subplot(1,2,2);
        dr = sqrt((p.dE*p.ev1x/2)^2+(p.dE*p.ev1y/2)^2);
        dC = dr/10;
        centers = (0:dC:(dr-dC)) + 0.5*dC;
        cE = hist(distE,centers);
        cI = hist(distI,centers);
        bar(centers,[cI', cE']);
        legend({'inh','exc'});
        xlabel('\mu m');
        ylabel('# of connection');

        if ~isempty(format)
            if strcmp(format,'fig')
                format = format0;
                printDriver = printDriver0;
            end
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(hEx,[outputfdr,'/',num2str(i),'-D-',matName,'.',format0],printDriver,dpi);
        end
    % I;
    else
        i = i-p.nv1e;
        subplot(1,2,1);
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
        title([num2str(i+p.nv1e),'th neuron, Inh, nLGN = ',num2str(nLGN(i+p.nv1e,1))]);
        xlim([min(min(pex),min(pix)),max(max(pex),max(pix))]);
        ylim([min(min(pey),min(piy)),max(max(pey),max(piy))]);   
        axis equal
    %   subplot(2,3,6);
        subplot(1,2,2);
        cE = hist(distE,centers);
        cI = hist(distI,centers);
        cR = hist(distE(retro),centers);
        bar(centers,[cI', cE',cR']);
        xlabel('\mu m');
        ylabel('# of connection');
        legend({'inh','exc','reciprocal'});
        if ~isempty(format)
            if strcmp(format,'fig')
                format = format0;
                printDriver = printDriver0;
            end
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(hEx,[outputfdr,'/',num2str(i+p.nv1e),'-D-',matName,'.',format0],printDriver,dpi);
        end
        end
    end
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
