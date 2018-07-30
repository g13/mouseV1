function plotFeedBackInh2(theme0,theme1,theme2,theme3,pfile,ntheta,format)
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
    contrastRange = [2,4];
    gIr = zeros(p.nv1,2,2);
    gI0 = zeros(p.nv1,2,2);
    gI1 = zeros(p.nv1,2,2);
    theme = {{theme0,theme1},{theme2,theme3}};
    eps = 125*2^(4-1);
    for ii = 1:2
        disp(['loading theme ',num2str(ii)]);
        DIR = [theme{ii}{1},'/',num2str(eps,'%04d')];
        [ipA0,ioA0,gi0] = readInh(DIR,ntheta,p.nv1);

        DIR = [theme{ii}{2},'/',num2str(eps,'%04d')];
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
    pick = 1:p.nv1e;
    [ax,h1,h2] = plotyy([1,2,3,4],[mean(gIr(:,1,1)),mean(gIr(:,2,1)),mean(gIr(:,1,2)),mean(gIr(:,2,2))],[1,2,3,4],zeros(4,1)-1);
    set(ax(1),'XTick',[1.5,3.5],'XTickLabel',{'0.6','\infty'});
    set(ax(1),'XLabel','\sigma(P_{I\to E})');
    h1.LineStyle = 'None';
    h1.Marker = '*';
    h1.MarkerSize = 8;
    h1.Color = 'r';
    xlim(ax(1),[0,5]);
    ylim(ax(1),[0,1]);
    set(ax(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
    ylabel(ax(1),'Feedback %');
    set(ax(1),'YColor','k');
    hold(ax(2))
    for i=1:2
        %errorbar(ax(2),2*(i-1)+1,mean(gI1(:,1,i)),std(gI1(:,1,i)),'b');
        bar(ax(2),2*(i-1)+1,mean(gI1(pick,1,i)),'b');
        %errorbar(ax(2),2*(i-1)+1,mean(gI0(:,1,i)),std(gI0(:,1,i)),'c');
        bar(ax(2),2*(i-1)+1,mean(gI0(pick,1,i)),'c');
        %errorbar(ax(2),2*(i-1)+2,mean(gI1(:,2,i)),std(gI1(:,2,i)),'b');
        bar(ax(2),2*(i-1)+2,mean(gI1(pick,2,i)),'b');
        %errorbar(ax(2),2*(i-1)+2,mean(gI0(:,2,i)),std(gI0(:,2,i)),'c');
        bar(ax(2),2*(i-1)+2,mean(gI0(pick,2,i)),'c');
    end
    set(ax(2),'YColor','b');
    set(ax(2),'XTick',[]);
    set(ax(1),'Color','none');
    set(ax(2),'Color','none');
    uistack(ax(1),'up',2);
    xlim(ax(2),[0,5]);
    ylim(ax(2),[0,inf]);
    set(ax(2), 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
    ylabel(ax(2),'g_{I}');

    if ~isempty(format)
        fname = ['FeedBackInh-',theme0,'-',theme1,'-',theme2,'-',theme3];
        saveas(h,fname);
        print(h,[fname,'.',format],printDriver,dpi);
    end
end
function [ipA,ioA,gi] = readInh(DIR,ntheta,n)
    filepath = [DIR, '/cv.dat'];
    fid = fopen(filepath);
    Data = fread(fid, [n, 16],'double');
    ipA = nearest(Data(:,15));
    ioA = nearest(Data(:,16));
    fclose(fid);

    filepath = [DIR, '/intra.dat'];
    fid = fopen(filepath);
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    fread(fid, [n, 2*ntheta+1],'double');
    gi = fread(fid, [n, 2*ntheta+1],'double');
    fclose(fid);
end
