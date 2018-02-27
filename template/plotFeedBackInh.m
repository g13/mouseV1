function plotFeedBackInh(theme0,theme1,lgnfile,ntheta,format)
    if nargin < 5
        format = '';
        if nargin < 4
            ntheta = 12;
        end
    end
    load(pfile,'p');
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r300';
    end
    dir = [dir,'/'];
    contrastRange = [2,4];
    gIr = zeros(p.nv1,2,length(contrastRange));
    gI0 = zeros(p.nv1,2,length(contrastRange));
    gI1 = zeros(p.nv1,2,length(contrastRange));
    for ii = 1:length(contrastRange)
        i = contrastRange(ii);
        disp(['loading c',num2str(i)]);
        eps = 125*2^(i-1);
        DIR = [theme0,'/',num2str(eps,'%04d')];
        [ipA0,ioA0,gi0] = readInh(DIR,ntheta,p.nv1);

        DIR = [theme1,'/',num2str(eps,'%04d')];
        [ipA1,ioA1,gi1] = readInh(DIR,ntheta,p.nv1);
        disp('loaded, recalibrating');
        for k=1:p.nv1
            gI0(k,1,ii) = gi0(k,ipA0(k));
            gI0(k,2,ii) = gi0(k,ioA0(k));

            gI1(k,1,ii) = gi1(k,ipA1(k));
            gI1(k,2,ii) = gi1(k,ioA1(k));

            gIr(k,1,ii) = (gi1(k,ipA0(k))-gi0(k,ipA0(k)))/gi1(k,ipA0(k));
            gIr(k,2,ii) = (gi1(k,ioA0(k))-gi0(k,ioA0(k)))/gi1(k,ioA0(k));
        end
        disp('recalibrated');
    end
    h = figure;
    subplot(2,1,1)
    hold on
    for i=1:length(contrastRange)
        errorbar(2*(i-1)+1,mean(gI0(:,1,i)),std(gI0(:,1,i)),'c');
        bar(2*(i-1)+1,mean(gI0(:,1,i)),'c');
        errorbar(2*(i-1)+1,mean(gI1(:,1,i)),std(gI1(:,1,i)),'b');
        bar(2*(i-1)+1,mean(gI1(:,1,i)),'b');
        errorbar(2*(i-1)+2,mean(gI0(:,2,i)),std(gI0(:,2,i)),'c');
        bar(2*(i-1)+2,mean(gI0(:,2,i)),'c');
        errorbar(2*(i-1)+2,mean(gI1(:,2,i)),std(gI1(:,2,i)),'b');
        bar(2*(i-1)+2,mean(gI1(:,2,i)),'b');
    end
    set(ax,'XTick',[1,2,3,4],'XTickLabel',{'o25%','p25%','o100%','p100%'});
    ax = subplot(2,1,2)
    plot([1,2,3,4],[mean(gIr(:,1,1)),mean(gIr(:,2,1)),mean(gIr(:,1,2)),mean(gIr(:,2,2))]);
    xlim([0,5]);
    set(ax,'XTick',[1,2,3,4],'XTickLabel',{'o25%','p25%','o100%','p100%'});

    if ~isempty(format)
        fname = [theme,'/FeedBackInh-',theme];
        saveas(h,fname);
        print(h,[fname,'.',format],printDriver,dpi);
    end
end
function [ipA,ioA,gi] = readInh(DIR,ntheta,n)
    filepath = [DIR, '/', 'intra.dat'];
    fid = fopen(filepath);
    Data = fread(fid, [n, 16],'double');
    ioA = nearest(Data(:,15));
    ipA = nearest(Data(:,16));

    filepath = [DIR, '/', 'intra.dat'];
    fid = fopen(filepath);
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    gi = fread(fid, [n, 2*ntheta+1],'double');
end
