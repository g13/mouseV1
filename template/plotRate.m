function plotRate(lgnfile,theme,eps,ntheta,format)
    pkEdgeE = [];
    pkEdgeI = [];
    rateEdgeE = [];
    rateEdgeI = [];
    spontEdgeE = [];
    spontEdgeI = [];
    if nargin < 5
        format = '';
        if nargin < 4
            ntheta = 12;
            if nargin < 3
                eps = 4;
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
        dpi = '-r300';
    end
    load(lgnfile,'p');
    epick = 1:p.nv1e;
    ipick = (p.nv1e+1):p.nv1;
    eps = 125*2^(eps-1);
    epsStr = num2str(eps,'%04d');
    DIR = [theme,'/',epsStr];
    [spontRate, pkRate, rate] = readRate(DIR,ntheta,p.nv1e,25);
    h = figure;
    subplot(3,2,1);
    if ~isempty(pkEdgeE)
        histogram(pkRate(epick),pkEdge);
    else
        histogram(pkRate(epick));
    end
    title('Excitatory');
    ylabel('# of Cells');
    xlabel('Peak Rate (Hz)');
    subplot(3,2,2);
    if ~isempty(pkEdgeI)
        histogram(pkRate(ipick),pkEdgeI);
    else
        histogram(pkRate(ipick));
    end
    title('Inhibitory');
    ylabel('# of Cells');
    xlabel('Peak Rate (Hz)');

    ax = subplot(3,2,3);
    barRate(rate(epick,1:ntheta),ax,rateEdgeE);
    ylabel('\bar{#} of Cells per orientation');
    xlabel('Rate (Hz)');
    ax = subplot(3,2,4);
    barRate(rate(ipick,1:ntheta),ax,rateEdgeI);
    ylabel('\bar{#} of Cells per orientation');
    xlabel('Rate (Hz)');

    subplot(3,2,5)
    histogram(spontRate(epick),spontEdgeE);
    ylabel('# of Cells');
    xlabel('Spont. Rate (Hz)');
    subplot(3,2,6)
    histogram(spontRate(ipick),spontEdgeI);
    ylabel('# of Cells');
    xlabel('Spont. Rate (Hz)');

    if ~isempty(format)
        fname = [theme,'/Rate-',epsStr,'-',theme];
        saveas(h,fname);
        print(h,[fname,'.',format],printDriver,dpi);
    end
end
function barRate(d,ax,edge,c,l,u)
    if nargin < 6
        u = 0.75;
        if nargin < 5
            l = 0.25;
            if nargin < 4
                c = [];
                if nargin < 3
                    edge = [];
                end
            end
        end
    end
    if ~isempty(edge)
        [N,~] = histcounts(d(:),edge);
    else
        [N,edge] = histcounts(d(:));
    end
    dedge = edge(2)-edge(1);
    bar(edge(1:length(edge)-1)+dedge/2,N/size(d,2));
end
