function plotFFT(pfile,theme,eps,dtheta,format)
    if nargin < 5
        format = '';
    end
    load(pfile,'p');
    pPosition = [0, 0, 1280, 720];
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r100';
    end
    n = p.nv1;
    FontSize = 20;
    set(0,'DefaultAxesFontSize',FontSize);
    
    fid = fopen([theme,'/sample_id.dat']);
    sampleID = fread(fid,[1,inf], 'int');
    nsample = length(sampleID);
    disp([num2str(nsample)]);
    fclose(fid);

    fid = fopen([theme,'/1000/cv.dat']);
    tmp = fread(fid,[n,16],'double');
    ipriA = tmp(:,15);
    fclose(fid);


    filepath = [theme,'/samples/1000/00-samples.dat'];
    fid = fopen(filepath)
    tmp = fread(fid,[9*nsample,inf],'double');
    fclose(fid);
    nt = size(tmp,2)
    clear tmp;
    Fs = 1/1e-3;
    v = zeros(nt,nsample,4);
    vs = zeros(nt,nsample,4);
    for ieps = 1:eps
        epsfdr = num2str(0.125*2^(ieps-1)*1000,'%04d');
        for is=1:nsample
            theta = ipriA(sampleID(is)) - dtheta;
            if theta < 0
                theta = 12 + theta;
            end
            if theta > 11
                theta = theta - 12;
            end
            thetafdr = num2str(theta,'%02d');
            
            filepath = [theme,'/samples/',epsfdr,'/',thetafdr,'-samples.dat'];
            fid = fopen(filepath);
            data = fread(fid,[9*nsample,inf],'double');
            fclose(fid);
            data = reshape(data,[9,nsample,nt]);

            v(:,is,ieps) = squeeze(data(1,is,:))';
            vs(:,is,ieps) = squeeze(data(8,is,:)./data(7,is,:))';
        end
    end
    h = figure; 
    epsRange = [2,4];
    for iieps = 1:2
        ieps = epsRange(iieps);
        ax = subplot(2,2,iieps);
        hold(ax);
        lineVE = FFTavg(v(:,sampleID<=p.nv1e,ieps),Fs,ax,200);
        lineVsE = FFTavg(vs(:,sampleID<=p.nv1e,ieps),Fs,ax,200);
        lineVE.Color = 'b';
        lineVsE.Color = 'r';
        legend({'V_{E}','Vs_{E}'});
        title(['C',num2str(12.5*2^(ieps-1)),'%']);
        xlabel('freq (Hz)');
        xlim([0,inf]);
        if ieps == 1
            ylabel('Exc Mem Spectrum');
        end
        
        ax = subplot(2,2,2+iieps);
        hold(ax);
        lineVI = FFTavg(v(:,sampleID>p.nv1e,ieps),Fs,ax,200);
        lineVsI = FFTavg(vs(:,sampleID>p.nv1e,ieps),Fs,ax,200);
        lineVI.Color = 'b';
        lineVsI.Color = 'r';
        legend({'V_{I}','Vs_{I}'});
        xlabel('freq (Hz)');
        xlim([0,inf]);
        if ieps == 1
            ylabel('Inh Mem Spectrum');
        end
    end
    if ~isempty(format)
        figname = [theme,'/d',num2str(dtheta),'-VsVsample'];
        saveas(h,figname);
        if ~strcmp(format,'fig')
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(h,[figname,'.',format],printDriver,dpi);
        end
    end
end
function line = FFTavg(d,Fs,ax,nf)
    n = size(d,2);
    nt = size(d,1);
    nfft = 2^nextpow2(nt);
    if nargin< 4
        nf = nfft/2+1;
    end
    Y = zeros(nf,n);
    f = Fs/2*linspace(0,1,nfft/2+1);
    f = f(1:nf);
    for i = 1:n
        y = fft(d(:,i)-mean(d(:,i)),nfft)/nt;
        Y(:,i) = 2*abs(y(1:nf));
    end
    line = plot(ax,f,mean(Y,2),'LineWidth',2);
    xlabel(ax,'Hz');
    xlim(ax,[0,inf]);
end
