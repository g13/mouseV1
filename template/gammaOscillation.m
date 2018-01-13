function gammaOscillation(theme,eps,theta,twindow,format)
    if nargin < 5
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
    
    FontSize = 20;
    set(0,'DefaultAxesFontSize',FontSize);
    fid = fopen([theme,'/theta_wise/',eps,'/',theta,'-gammaOsci.dat']);
    data = fread(fid,[2,inf],'float');
    fclose(fid);
    l = size(data,2);
    dt = twindow/(l+1);
    t = linspace(dt,twindow-dt,l)*1000;
    
    Fs = 1/1e-4;
    NFFT = 2^nextpow2(l);
    sp = data(1,:);
    psp = data(2,:);
    
    h = figure; 
    subplot(2,2,1);
    ssp = smooth(sp,81);
    plot(t,ssp);
    title([theme,'-',eps,'-',theta]);
    ylabel('summed spikes');
    xlabel('t (ms)');
    xlim([0,200]);
    subplot(2,2,2);
    hold on
    grid on
    Y = fft(ssp-mean(ssp),NFFT)/l;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    Y = 2*abs(Y(1:NFFT/2+1));
    plot(f,Y);
    [m,i] = max(Y);
	plot(f(i),m,'*r');
  %  sY = smooth(Y,11);
    
   % plot(f,sY,'-r','LineWidth',3);
    xlabel('Hz');
    title('FFT');
    xlim([0,500]);
    
    subplot(2,2,3);
    plot(t,psp);
    ylabel('agregated PSP');
    xlabel('t /ms');
    
    subplot(2,2,4);
    hold on
    grid on
    Y = fft(psp-mean(psp),NFFT)/l;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    Y = 2*abs(Y(1:NFFT/2+1));
    plot(f,Y);
    [m,i] = max(Y);
	plot(f(i),m,'*r');
    xlabel('Hz');
    xlim([0,500]);
    if ~isempty(format)
        figname = [theme,'/',eps,'-',theta,'-gammaOsci.',format];
        if strcmp(format,'fig')
            saveas(h,figname);
        else
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(h,figname,printDriver,dpi);
        end
    end
end
