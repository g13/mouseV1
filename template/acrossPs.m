function acrossPs(lgnfile,x,fdr,ntheta,outputfdr,theme,format)
    nfdr = length(fdr);
    assert(length(x)==nfdr);
    if nargin < 7
        format = '';
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
    FontSize = 10;
    legendFontSizeOffset = 2;
    set(0,'DefaultAxesFontSize',FontSize)
%dataRead = true;
dataRead = false;
% readData
load(lgnfile);
nfdr = length(fdr)
cl = 4;
sc100 = zeros(nfdr,p.nv1);
cv100 = sc100;
cv50 = sc100;
for i = 1:nfdr
    load([fdr{i},'-tcData-x',num2str(cl),'.mat'],'nP');
    sc100(i,:) = nP(4).sc;
    cv50(i,:) = nP(2).cv;
    cv100(i,:) = nP(4).cv;
    clear nP 
end
%% CV % above y=x line
hCVabove = figure;
cvabe = zeros(nfdr,1);
cvbei = zeros(nfdr,1);

cvave = zeros(nfdr,1);
cvavi = zeros(nfdr,1);
dcvave = zeros(nfdr,2);
dcvavi = zeros(nfdr,2);

for i = 1:nfdr
    temp =  cv100(i,:)>cv50(i,:);
    cvabe(i) = sum(temp(1:p.nv1e))/p.nv1e*100;
    cvbei(i) = sum(~temp(p.nv1e+(1:p.nv1i)))/p.nv1i*100;
    
    temp = cv100(i,:)-cv50(i,:);
    dcvave(i,1)  = mean(temp(1:p.nv1e));
    dcvave(i,2)  = std(temp(1:p.nv1e));
    dcvavi(i,1)  = mean(-temp(p.nv1e+(1:p.nv1i)));
    dcvavi(i,2)  = std(temp(p.nv1e+(1:p.nv1i)));
    
    cvave(i) = mean(cv100(i,1:p.nv1e));
    cvavi(i) = mean(cv50(i,p.nv1e+(1:p.nv1i)));
end
clear temp
subplot(1,2,1);
hold on
plot(x,cvabe,'or');
plot(x,cvbei,'ob');
legend({'e','i'});
title('CV over y=x')
ylabel('% of Population')
title(theme);
subplot(2,2,2);
bar(x,dcvave);
hold on
plot(x,cvave,'or');

%legend({'\Delta CV(25% vs 100%)','absolute CV(100%)'});
title('Exc');
subplot(2,2,4);
bar(x,dcvavi);
hold on
plot(x,cvavi,'ob');
legend({'\Delta CV(25% vs 100%)','absolute CV(25%)'});
title('Inh');

if ~isempty(format)
    set(gcf,'Renderer','Painters')
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
    print(hCVabove,[outputfdr,'/',theme,'-','-pP-CVabove.',format],printDriver,dpi);
end

%% F1/F0 < 1 % in ORF and Inh peak positions 

hSC = figure;
peakPos = zeros(i,1)
per
for



if ~isempty(format)
    set(gcf,'Renderer','Painters')
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
    print(hSC,[outputfdr,'/',theme,'-','-pP-SC.',format],printDriver,dpi);
end
end
