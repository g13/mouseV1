function FFTsample(theme,eps,dtheta,format)
    if nargin < 4
        format = '';
    end
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
    n = 10800; 
    FontSize = 20;
    set(0,'DefaultAxesFontSize',FontSize);
    
    fid = fopen([theme,'/sample_id.dat']);
    sampleID = fread(fid,[1,inf], 'int');
    nsample = length(sampleID);
    disp([num2str(nsample)]);
    fclose(fid);

    epsfdr = num2str(0.125*2^(eps-1)*1000,'%04d');
    fid = fopen([theme,'/',epsfdr,'/cv.dat']);
    tmp = fread(fid,[n,16],'double');
    ipriA = tmp(:,15);
    fclose(fid);


    filepath = [theme,'/samples/',epsfdr,'/00-samples.dat'];
    fid = fopen(filepath)
    tmp = fread(fid,[9*nsample,inf],'double');
    fclose(fid);
    nt = size(tmp,2)
    clear tmp;
    t = (0:(nt-1));
    NFFT = 2^nextpow2(nt);
    Fs = 1/1e-3;
    sumVfft = zeros(NFFT/2+1,1);
    sumVsfft = zeros(NFFT/2+1,1);
    h = figure; 
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
        fid = fopen(filepath)
        data = fread(fid,[9*nsample,inf],'double');
        fclose(fid);
        data = reshape(data,[9,nsample,nt]);

        vi = squeeze(data(1,is,:));
        vsi = squeeze(data(8,is,:)./data(7,is,:));
        
        subplot(2,2,1);
        hold on
        %ssp = smooth(v,81);
        plot(t,vi);
        plot(t(nt),mean(vi),'*r');
        %title([theme,'-c',num2str(eps),'_d',num2str(dtheta)]);
        if is == nsample
            ylabel('Mem');
            xlabel('t (ms)');
            xlim([0,1000]);
        end
        subplot(2,2,2);
        hold on
        grid on
        Y = fft(vi-mean(vi),NFFT)/nt;
        f = Fs/2*linspace(0,1,NFFT/2+1);
        Y = 2*abs(Y(1:NFFT/2+1));
        sumVfft = sumVfft + Y;
        plot(f,Y);
  %      sY = smooth(Y,11);
        
   %     plot(f,sY,'-r','LineWidth',3);
        if is == nsample
            plot(f,sumVfft/nsample,'k','LineWidth',2);
            xlabel('Hz');
            title('Mem FFT');
            xlim([0,inf]);
        end
        
        subplot(2,2,3);
        hold on
        plot(t,vsi);
        plot(t(nt),mean(vsi),'*r');
        if is == nsample
            xlim([0,1000]);
            ylabel('Vs');
            xlabel('t (ms)');
        end
        subplot(2,2,4);
        hold on
        grid on
        Y = fft(vsi-mean(vsi),NFFT)/nt;
        f = Fs/2*linspace(0,1,NFFT/2+1);
        Y = 2*abs(Y(1:NFFT/2+1));
        sumVsfft = sumVsfft + Y;
        plot(f,Y);
        if is == nsample
            plot(f,sumVsfft/nsample,'k','LineWidth',2);
            title('Vs FFT');
            xlabel('Hz');
            xlim([0,inf]);
        end
    end
    if ~isempty(format)
        figname = [theme,'/c',num2str(eps),'_d',num2str(dtheta),'-VsVsample.',format];
        saveas(h,figname);
        if ~strcmp(format,'fig')
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(h,figname,printDriver,dpi);
        end
    end
end
