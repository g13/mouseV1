function plotresult(theme,lgn,ntheta,contrastLevel,format,outputfdr,thres)
addpath(genpath('../matlab_Utilities'));
if nargin < 7
	thres = 0.3;
	if nargin < 6
	    outputfdr = theme;
	    if nargin < 5
	        format = '';
	        if nargin < 4
	            contrastLevel = 4;
	            if nargin < 3
	                ntheta = 12;
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
		dpi = '-r100';
	end
LineWidth = 2;
set(groot,'defaultLineLineWidth',2);
FontSize = 14;
set(groot,'defaultAxesFontSize',FontSize);
set(groot,'defaultTextFontSize',FontSize);
set(groot,'defaultLegendFontSize',FontSize-1);
% position =[960,150,900,600];
pPosition = [0,0, 1280, 720];
position = [50,50,1400,900];
set(0, 'DefaultFigurePosition', position);
% theme = 'C7';
% theme = 'rv1';
load(lgn,'p','nLGN');
nLGN = nLGN(:,1);
for ifdr = 1:contrastLevel
eps = 125*2^(ifdr-1);
% eps = 250*ifdr;
fdr = [theme,'/',num2str(eps,'%0.4d')];
% disp(['contrast = ',num2str(eps/10,'%3.1f'),'%']);
% fdr = 'nm-04';
% fdr = 'mpin-02';
pw = false;
% pw = true;
% plotNoBack = true;
plotNoBack = false;

careBrate = true;
f4=fopen([fdr,'/cv.dat'],'r');

cv1d=fread(f4,p.nv1,'double');
cv1d = 1-cv1d;
sc1d=fread(f4,p.nv1,'double');
% sc1d = sc1d*2;
sc2_1d=fread(f4,p.nv1,'double');
mrate1d=fread(f4,p.nv1,'double');
pkrate1d=fread(f4,p.nv1,'double');
porate1d=fread(f4,p.nv1,'double');
po1d=fread(f4,p.nv1,'double');
excite1d=fread(f4,p.nv1,'double');
brate1d=fread(f4,p.nv1,'double');
cvNoBack1d=fread(f4,p.nv1,'double');
fclose(f4);

lgnmax_e = max(nLGN(excite1d>0.5));
lgnmin_e = min(nLGN(excite1d>0.5));
lgnmax_i = max(nLGN(excite1d<0.5));
lgnmin_i = min(nLGN(excite1d<0.5));

% f1 = fopen([fdr,'/intraavg.dat'],'r');
% fread(f1,[17,p.nv1],'double');
% vslaveSC = fread(f1,p.nv1,'double');
% fclose(f1);

%f1=fopen([fdr, '/08/i-and-f.dat3'],'r');
%see=fread(f1,p.nv1,'double');
%sei=fread(f1,p.nv1,'double');
%sie=fread(f1,p.nv1,'double');
%sii=fread(f1,p.nv1,'double');
%excite1d=fread(f1,p.nv1,'double');
%nlgn1d=fread(f1,p.nv1,'double');
%fexc1d=fread(f1,p.nv1,'double');
%finh1d=fread(f1,p.nv1,'double');
%fclose(f1);
%ce0=6.0;ci0=150;
%feee=fexc1d*ce0;
%fiii=finh1d*ci0;
%feee0=0.5d0*feee;
%fiii0=0.5d0*fiii;


if careBrate
    efired=excite1d > .5 & pkrate1d > brate1d & pkrate1d > thres;
    ifired=excite1d < .5 & pkrate1d > brate1d & pkrate1d > thres;
else
    efired=excite1d > .5 & pkrate1d > thres;
    ifired=excite1d < .5 & pkrate1d > thres;
end
range = 2.0/20.:2.0/10.:2.0-2.0/20.;
edges = 0:0.2:2;

nc=0; ns=0;
bratec=[];brates=[];

for i=1:p.nv1/2,
    if sc1d(i) > 1
        ns = ns + 1; 
        brates(ns) = brate1d(i);
    else
        nc = nc + 1; 
        bratec(nc)=brate1d(i);
    end
end


% figure;
% subplot(2,1,1);
% y=hist(vslaveSC(efired),range);bar(0.05:0.2:0.95,y(1:5),'b');hold on;
% bar(1.05:0.2:1.95,y(6:10),'r');hold off;
% set(gca,'FontSize',FontSize);set(gca,'FontName','Times');
% % title('Histogram of V.F1/F0 ');
% ylabel('No. of Exc Cells');xlabel('V.F1/F0');
% xlim([0 2]);
% ylim([0,inf]);
% 
% subplot(2,1,2);
% y=hist(vslaveSC(ifired),range);bar(0.05:0.2:0.95,y(1:5),'b');hold on;
% bar(1.05:0.2:1.95,y(6:10),'r');hold off;
% set(gca,'FontSize',FontSize);set(gca,'FontName','Times');
% % title('Histogram of V.F1/F0 ');
% ylabel('No. of Inh Cells');xlabel('V.F1/F0');
% xlim([0 2]);
% ylim([0,inf]);

%efired=find(excite1d > .5 & mrate1d > 8);
%ifired=find(excite1d < .5 & mrate1d > 8);
h = figure;
subplot(4,3,1);
%y=hist(sc1d(efired),range);
y=histcounts(sc1d(efired),edges);
bar(0.1:0.2:0.9,y(1:5),'b');hold on;
bar(1.1:0.2:1.9,y(6:10),'r');hold off;
title('Histogram of F1/F0');ylabel('# Exc');
set(gca,'FontSize',FontSize);
% xlabel('F1/F0');
xlim([0 2]);
ylim([0,inf]);

subplot(4,3,4);
y=hist(sc1d(ifired),range);
bar(0.1:0.2:0.9,y(1:5),'b');hold on;
bar(1.1:0.2:1.9,y(6:10),'r');hold off;
ylabel('# Inh');xlabel('F1/F0');
set(gca,'FontSize',FontSize);
xlim([0 2]);
ylim([0,inf]);

subplot(6,3,2);y=hist(cv1d(ifired),0.05:0.1:0.95);bar(0.05:0.1:0.95,y,'k');
axis([0 1 0 inf]);
set(gca,'FontSize',FontSize);
set(gca,'XTickLabel',['']);
title('Histogram of 1-CV, Inh');

if careBrate
    %ec=find(excite1d  > .5 & sc1d < .5 & mrate1d > 8);
    %es=find(excite1d > .5 & sc1d > .5 & mrate1d > 8);
    ec=find(excite1d  > .5 & sc1d <1.0 & pkrate1d > brate1d & pkrate1d > thres);
    es=find(excite1d > .5 & sc1d > 1.0 & pkrate1d > brate1d & pkrate1d > thres);
else
    ec=find(excite1d > .5 & sc1d < 1.0 & pkrate1d > thres);
    es=find(excite1d > .5 & sc1d > 1.0 & pkrate1d > thres);
end
cvs=cv1d(es);
cvc=cv1d(ec);
brates=brate1d(es);
bratec=brate1d(ec);

subplot(6,3,5);y=hist(cvs,0.05:0.1:0.95);
bar(0.05:0.1:0.95,y,'r'); axis([0 1 0 inf]);
set(gca,'FontSize',FontSize);
set(gca,'XTickLabel',['']);
ylabel('# of Cells');
title('Simple Exc');
% legend('Exc. Simple',2);

subplot(6,3,8);
y=hist(cvc,0.05:0.1:0.95);
bar(0.05:0.1:0.95,y,'b'); axis([0 1 0 inf]);
title('Complex Exc');
xlabel('1-CV');
set(gca,'FontSize',FontSize);
% legend('Exc. Complex',2);

% subplot(6,3,3);
% hist(brate1d(ifired));
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','k');
% % plot(cv1d(ifired),brate1d(ifired),'ko');
% set(gca,'FontSize',FontSize);
% title('Histogram of Spontaneous Rate');
% 
% % title('CV vs. Spontaneous Rate');
% legend('Inh',1);
% % % axis([0 1 0 5]);
% axis tight
% 
% subplot(6,3,6);
% hist(brates);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r');
% % plot(cvs,brates,'ro');
% set(gca,'FontSize',FontSize);
% ylabel('# Neurons');
% % ylabel('Spontaneous Firing Rate (sec^{-1})');
% legend('Exc. Simple',1);
% % % axis([0 1 0 5]);
% axis tight
% 
% subplot(6,3,9);
% hist(bratec);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b');
% % plot(cvc,bratec,'bo');
% set(gca,'FontSize',FontSize);
% legend('Exc. Complex',1);
% % axis([0 1 0 5]);
% xlabel('Spontaneous Rate (Hz)');
% axis tight


subplot(4,3,3);hold on;

ylabel('Spont. Rate (Hz)');
errorbar(1:4,[mean(brate1d(excite1d>.5)),mean(brate1d(efired)),mean(brate1d(excite1d<.5)),mean(brate1d(ifired))],...
            [std(brate1d(excite1d>.5)),std(brate1d(efired)),std(brate1d(excite1d<.5)),std(brate1d(ifired))],'.');
bar([mean(brate1d(excite1d>.5)),mean(brate1d(efired)),mean(brate1d(excite1d<.5)),mean(brate1d(ifired))]);
set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'','','',''});

subplot(4,3,6);hold on
ylabel('Evoked. Rate (Hz)');
errorbar(1:4,[mean(pkrate1d(excite1d>.5)),mean(pkrate1d(efired)),mean(pkrate1d(excite1d<.5)),mean(pkrate1d(ifired))],...
            [std(pkrate1d(excite1d>.5)),std(pkrate1d(efired)),std(pkrate1d(excite1d<.5)),std(pkrate1d(ifired))],'.');
bar([mean(pkrate1d(excite1d>.5)),mean(pkrate1d(efired)),mean(pkrate1d(excite1d<.5)),mean(pkrate1d(ifired))]);
set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'Exc','E-active','Inh','I-active'});

subplot(4,2,5);
hist(po1d(efired)/pi*180,ntheta);
xlabel('Exc preferred angle');
ylabel('# Exc');
xlim([0,180-180/ntheta]);

subplot(4,2,7);
hist(po1d(ifired)/pi*180,ntheta);
xlabel('Inh preferred angle');
ylabel('# Inh');
xlim([0,180-180/ntheta]);

subplot(2,2,4);
hold on;
for j=min(lgnmin_e,lgnmin_i):max(lgnmax_e,lgnmax_i)
    if careBrate
        pick0 = excite1d > .5 & pkrate1d > brate1d & nLGN==j & pkrate1d > thres;
    else
        pick0 = excite1d > .5 & nLGN==j & pkrate1d > thres;
    end
    if ~isempty(pick0)
        pick = pick0 & [p.typeE==1;false(p.nv1i,1)];
        if ~isempty(pick)
            if ~plotNoBack
                plot3(cv1d(pick),sc1d(pick),pkrate1d(pick),'o','Color',[0.1+(j-lgnmin_e)/(lgnmax_e-lgnmin_e)*0.9,0,0]);
            else
                plot3(cvNoBack1d(pick),sc1d(pick),pkrate1d(pick),'o','Color',[0.1+(j-lgnmin_e)/(lgnmax_e-lgnmin_e)*0.9,0,0]);
            end
        end
        pick = pick0 & [p.typeE==2;false(p.nv1i,1)];
        if ~isempty(pick)
            if ~plotNoBack
                plot3(cv1d(pick),sc1d(pick),pkrate1d(pick),'o','Color',[0,0.1+(j-lgnmin_e)/(lgnmax_e-lgnmin_e)*0.9,0]);
            else
                plot3(cvNoBack1d(pick),sc1d(pick),pkrate1d(pick),'o','Color',[0,0.1+(j-lgnmin_e)/(lgnmax_e-lgnmin_e)*0.9,0]);
            end
        end
    end
    if careBrate
        pick = excite1d < .5 & pkrate1d > brate1d & nLGN==j & pkrate1d > thres;
    else
        pick = excite1d < .5 & nLGN==j & pkrate1d > thres;
    end                                
    if ~isempty(pick)
        if ~plotNoBack
            plot3(cv1d(pick),sc1d(pick),pkrate1d(pick),'*','Color',[0,0,0.1+(j-lgnmin_i)/(lgnmax_i-lgnmin_i)*0.9]);
        else
            plot3(cvNoBack1d(pick),sc1d(pick),pkrate1d(pick),'*','Color',[0,0,0.1+(j-lgnmin_i)/(lgnmax_i-lgnmin_i)*0.9]);
        end
    end
end
axis([0,1,0,2]);
campos([-3,-5,max(pkrate1d)*10]);
grid on;

set(gca,'FontSize',FontSize);
xlabel('CV');ylabel('F1/F0');zlabel('peak Evoked Rate');
title({fdr,['evoked E~',num2str(sum(efired/p.nv1e)*100,'%1.1f'),'% I~',num2str(sum(ifired/p.nv1i)*100,'%1.1f'),'%']});
if ~isempty(format)
	set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);	
    if strcmp(format,'fig')
        saveas(h,[outputfdr,'/',theme,'-',num2str(eps,'%0.4d'),'-population.',format]);
    else
        print(h,[outputfdr,'/',theme,'-',num2str(eps,'%0.4d'),'-population.',format],printDriver,dpi);
    end
end

hTypes = figure;
evoked = zeros(p.ntypeE,1);
for i=1:p.ntypeE
	evoked(i) = 100 * sum([p.typeE == i;false(p.nv1i,1)] & efired) / sum(p.typeE == i);	
end
evokedP = strcat('evoked:',cellstr(num2str(evoked,'%3.1f')),'%');

for i = 1:p.ntypeE
    pick = efired & [p.typeE == i; false(p.nv1i,1)];

    subplot(4, p.ntypeE + p.ntypeI, i)
    %y=hist(sc1d(pick),range);
    %bar(0.05:0.2:0.95,y(1:5),'b');hold on;
    %bar(1.05:0.2:1.95,y(6:10),'r');hold off;
    y=histcounts(sc1d(pick),edges);
    bar(0.1:0.2:0.9,y(1:5),'b');hold on;
    bar(1.1:0.2:1.9,y(6:10),'r');hold off;
    title(p.Etypes{i});

	xlabel('F1/F0'); ylabel('# Cell');
	xlim([0 2]);
    
    subplot(4, p.ntypeE + p.ntypeI, p.ntypeE + p.ntypeI + i)
	y=hist(cv1d(pick),0.05:0.1:0.95);
	bar(0.05:0.1:0.95,y,'r');
	xlabel('1-CV'); ylabel('# Cell');
	title(['mean ',num2str(mean(cv1d(pick)),'%1.2f')]);

    subplot(4, p.ntypeE + p.ntypeI, 2*p.ntypeE + 2*p.ntypeI + i)
	hist(pkrate1d(pick),10);
	%bar(0.05:0.1:0.95,y,'r');
	xlabel('FR'); ylabel('# Cell');
	title(['mean ',num2str(mean(pkrate1d(pick)),'%3.1f'),'Hz']);

    subplot(4, p.ntypeE + p.ntypeI, 3*p.ntypeE + 3*p.ntypeI + i)
	hist(po1d(pick)/pi*180,ntheta);
	xlabel('pref angle'); ylabel('# Cell');
	title(evokedP{i});
	xlim([0,180-180/ntheta]);
end

evoked = zeros(p.ntypeI,1);
for i=1:p.ntypeI
	evoked(i) = 100 * sum([false(p.nv1e,1);p.typeI == i] & ifired) / sum(p.typeI == i);	
end
evokedP = strcat('evoked:',cellstr(num2str(evoked,'%3.1f')),'%');

for i = 1:p.ntypeI
    pick = ifired & [false(p.nv1e,1); p.typeI == i];

    subplot(4, p.ntypeE + p.ntypeI, p.ntypeE + i)
	hold on;
    y=histcounts(sc1d(pick),edges);
    bar(0.1:0.2:0.9,y(1:5),'b');hold on;
    bar(1.1:0.2:1.9,y(6:10),'r');hold off;
    %y=hist(sc1d(pick),range);
    %bar(0.05:0.2:0.95,y(1:5),'b');
    %bar(1.05:0.2:1.95,y(6:10),'r');
    title(p.Itypes{i});
	xlabel('F1/F0'); ylabel('# Cell');
	xlim([0 2]);

    subplot(4, p.ntypeE + p.ntypeI, 2*p.ntypeE + p.ntypeI + i)
	y=hist(cv1d(pick),0.05:0.1:0.95);
	bar(0.05:0.1:0.95,y,'b');
	xlabel('1-CV'); ylabel('# Cell');
	title(['mean ',num2str(mean(cv1d(pick)),'%1.2f')]);
    
    subplot(4, p.ntypeE + p.ntypeI, 3*p.ntypeE + 2*p.ntypeI + i)
	hist(pkrate1d(pick),10);
	%bar(0.05:0.1:0.95,y,'b');
	xlabel('FR'); ylabel('# Cell');
	title(['mean ',num2str(mean(pkrate1d(pick)),'%3.1f'),'Hz']);

    subplot(4, p.ntypeE + p.ntypeI, 4*p.ntypeE + 3*p.ntypeI + i)
	hist(po1d(pick)/pi*180,ntheta);
	xlabel('pref angle'); ylabel('# Cell');
	xlim([0,180-180/ntheta]);
	title(evokedP{i});
end

if ~isempty(format)
	set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);	
    if strcmp(format,'fig')
        saveas(hTypes,[outputfdr,'/',theme,'-',num2str(eps,'%0.4d'),'-type_specific.',format]);
    else
        print(hTypes,[outputfdr,'/',theme,'-',num2str(eps,'%0.4d'),'-type_specific.',format],printDriver,dpi);
    end
end
if ifdr==contrastLevel
hOverlap = figure;
ctrs = cell(2,1);
nX = 10;
sX = 5;
nY0 = 10;
sY = 5;
maxDis = max(p.enormDistance);
minDis = min(p.enormDistance);
[tickX, nX] = autoAxis(minDis,maxDis,nX,[0,inf]);
dTickX = tickX(2)-tickX(1);
tickLabelX = num2str(tickX');
lctrsx = (nX-1)*sX+1;
ctrs{1} = linspace(tickX(1),tickX(end),lctrsx);
tickPosX = linspace(0,lctrsx-1,nX)+0.5;

h1 = subplot(2,2,1);
maxSC = max(sc1d);
minSC = min(sc1d);
sc1de = sc1d(excite1d>0.5);
if minSC ~= maxSC
    [tickY, nY] = autoAxis(minSC,maxSC,nY0,[0,2]);
    dTickY = tickY(2)-tickY(1);
    tickLabelY = num2str(flipud(tickY'));
    ctrs{2} = linspace(tickY(1),tickY(end),(nY-1)*sY+1);
    lctrsy = length(ctrs{2});
    tickPosY = linspace(0,lctrsy-1,nY)+0.5;
    
    pair = [p.enormDistance, sc1de];
    den = hist3(pair,'Edges',ctrs);
    den = den(1:lctrsx-1,1:lctrsy-1);
    mdeny = max(den);
    [~, mden] = max(mdeny);
    den = den/mden;
    imagesc([1,lctrsx-1],[lctrsy-1,1],den');
    set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
    xlabel('norm'' Distance');
    ylabel('F1/F0');
    colormap(h1,redOnly);

    subplot(2,2,3)
    [~,ind] = histc(p.enormDistance,ctrs{1});
    enD = zeros(lctrsx,2);
    for i = 1:lctrsx
        enD(i,1) = mean(sc1de(ind==i));
        enD(i,2) = std(sc1de(ind==i));
    end
    errorbar(ctrs{1},enD(:,1),enD(:,2));
    xlim([ctrs{1}(1),ctrs{1}(end)]);
    xlabel('norm'' Distance');
    ylabel('F1/F0');
else
    disp(['all CV = ', num2str(minSC)]);
end

h2 = subplot(2,2,2);
maxCV = max(cv1d);
minCV = min(cv1d);
cv1de = cv1d(excite1d>.5);
if minCV ~= maxCV
    [tickY, nY] = autoAxis(minCV,maxCV,nY0,[0,1]);
    dTickY = tickY(2)-tickY(1);
    tickLabelY = num2str(flipud(tickY'));
    ctrs{2} = linspace(tickY(1),tickY(end),(nY-1)*sY+1);
    lctrsy = length(ctrs{2});
    tickPosY = linspace(0,lctrsy-1,nY)+0.5;

    pair = [p.enormDistance, cv1de];
    den = hist3(pair,'Edges',ctrs);
    den = den(1:lctrsx-1,1:lctrsy-1);
    mdeny = max(den);
    [~, mden] = max(mdeny);
    den = den/mden;
    imagesc([1,lctrsx-1],[lctrsy-1,1],den');
    set(gca,'YTickLabel',tickLabelY,'YTick',tickPosY,'XTickLabel',tickLabelX,'XTick',tickPosX);
    xlabel('norm'' Distance');
    ylabel('CV');
    colormap(h2,redOnly);

    subplot(2,2,4)
    enD = zeros(lctrsx,2);
    [~,ind] = histc(p.enormDistance,ctrs{1});
    for i = 1:lctrsx
        enD(i,1) = mean(cv1de(ind==i));
        enD(i,2) = std(cv1de(ind==i));
    end
    errorbar(ctrs{1},enD(:,1),enD(:,2));
    xlim([ctrs{1}(1),ctrs{1}(end)]);
    xlabel('norm'' Distance');
    ylabel('CV');
else
    disp(['all CV = ', num2str(minSC)]);
end

if ~isempty(format)
	set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);	
    if strcmp(format,'fig')
        saveas(hOverlap,[outputfdr,'/',theme,'-',num2str(eps,'%0.4d'),'-Overlap.',format]);
    else
        print(hOverlap,[outputfdr,'/',theme,'-',num2str(eps,'%0.4d'),'-Overlap.',format],printDriver,dpi);
    end
end
end
%%
if pw
    fii=fopen([fdr,'/intra.dat'],'r');
    cond=fread(fii,[p.nv1*17*9],'double');fclose(fii);
    cond=reshape(cond,[p.nv1 17 9]);

    gt=cond(1:p.nv1,1:16,1);
    cr=cond(1:p.nv1,1:16,2);
    vs=cond(1:p.nv1,1:16,3);
    gl=cond(1:p.nv1,1:16,4);
    ge=cond(1:p.nv1,1:16,5);
    gi=cond(1:p.nv1,1:16,6);
    fe=cond(1:p.nv1,1:16,7);
    fi=cond(1:p.nv1,1:16,8);
    rate=cond(1:p.nv1,1:16,9);

    gtb=cond(1:p.nv1,17,1);
    crb=cond(1:p.nv1,17,2);
    vsb=cond(1:p.nv1,17,3);
    glb=cond(1:p.nv1,17,4);
    geb=cond(1:p.nv1,17,5);
    gib=cond(1:p.nv1,17,6);
    feb=cond(1:p.nv1,17,7);
    fib=cond(1:p.nv1,17,8);
    rateb=cond(1:p.nv1,1:7,9);

    fii=fopen([fdr,'/intraavg.dat'],'r');
    dd=fread(fii,[p.nv1 20],'double');fclose(fii);
    ginh=dd(1:p.nv1,6);
    gexc=dd(1:p.nv1,5);
    gtot=dd(1:p.nv1,1);
    glgn=dd(1:p.nv1,4);
    curr=dd(1:p.nv1,2);
    vslv=dd(1:p.nv1,3);
    gett=gexc+glgn;
    finh=dd(1:p.nv1,8);
    fexc=dd(1:p.nv1,7);
    ginh0=dd(1:p.nv1,14);
    gexc0=dd(1:p.nv1,13);
    gtot0=dd(1:p.nv1,9);
    glgn0=dd(1:p.nv1,12);
    curr0=dd(1:p.nv1,10);
    vslv0=dd(1:p.nv1,11);
    gett0=gexc0+glgn0;
    finh0=dd(1:p.nv1,16);
    fexc0=dd(1:p.nv1,15);
    vssc1d=dd(1:p.nv1,17);
    vssc2_1d=dd(1:p.nv1,18);
    vscv1d=dd(1:p.nv1,19);
    vscv1_1d=dd(1:p.nv1,20);

    vs1d=po1d; poind1d=(po1d+pi/8)/pi*8;
    for i=1:p.nv1,vs1d(i) = find(vs(i,1:8) == max(vs(i,:)),1);end
    pvs1d=(vs1d-1)*pi/8.;
    fii=fopen([fdr,'/intra2.dat'],'r');
    cond2=fread(fii,p.nv1*17*9,'double');
    cond2=reshape(cond2,[p.nv1 17 9]);fclose(fii);

    gtb2=cond2(1:p.nv1,17,1);
    crb2=cond2(1:p.nv1,17,2);
    vsb2=cond2(1:p.nv1,17,3);
    glb2=cond2(1:p.nv1,17,4);
    geb2=cond2(1:p.nv1,17,5);
    gib2=cond2(1:p.nv1,17,6);
    feb2=cond2(1:p.nv1,17,7);
    fib2=cond2(1:p.nv1,17,8);

    gt2=cond2(1:p.nv1,1:16,1);
    cr2=cond2(1:p.nv1,1:16,2);
    vs2=cond2(1:p.nv1,1:16,3);
    gl2=cond2(1:p.nv1,1:16,4);
    ge2=cond2(1:p.nv1,1:16,5);
    gi2=cond2(1:p.nv1,1:16,6);
    fe2=cond2(1:p.nv1,1:16,7);
    fi2=cond2(1:p.nv1,1:16,8);

    gtstdev=sqrt(abs(gt2 - gt.*gt));
    crstdev=sqrt(abs(cr2 - cr.*cr));
    vsstdev=sqrt(abs(vs2 - vs.*vs));
    glstdev=sqrt(abs(gl2 - gl.*gl));
    gestdev=sqrt(abs(ge2 - ge.*ge));
    gistdev=sqrt(abs(gi2 - gi.*gi));
    festdev=sqrt(abs(fe2 - fe.*fe));
    fistdev=sqrt(abs(fi2 - fi.*fi));
    %rate=cond(1:p.nv1,1:16,7)/8.;

    gtbstdev=sqrt(abs(gtb2 - gtb.*gtb));
    crbstdev=sqrt(abs(crb2 - crb.*crb));
    vsbstdev=sqrt(abs(vsb2 - vsb.*vsb));
    glbstdev=sqrt(abs(glb2 - glb.*glb));
    gebstdev=sqrt(abs(geb2 - geb.*geb));
    gibstdev=sqrt(abs(gib2 - gib.*gib));
    febstdev=sqrt(abs(feb2 - feb.*feb));
    fibstdev=sqrt(abs(fib2 - fib.*fib));

    fii=fopen([fdr,'/intraavg2.dat'],'r');
    dd=fread(fii,[p.nv1 16],'double');fclose(fii);
    ginh2=dd(1:p.nv1,6);
    gexc2=dd(1:p.nv1,5);
    gtot2=dd(1:p.nv1,1);
    glgn2=dd(1:p.nv1,4);
    curr2=dd(1:p.nv1,2);
    vslv2=dd(1:p.nv1,3);
    finh2=dd(1:p.nv1,8);
    fexc2=dd(1:p.nv1,7);
    gett2=gexc2+glgn2;
    ginh02=dd(1:p.nv1,14);
    gexc02=dd(1:p.nv1,13);
    gtot02=dd(1:p.nv1,9);
    glgn02=dd(1:p.nv1,12);
    curr02=dd(1:p.nv1,10);
    vslv02=dd(1:p.nv1,11);
    finh02=dd(1:p.nv1,16);
    fexc02=dd(1:p.nv1,15);
    gett02=gexc0+glgn0;
    distratecorr;
end
end
