function VsFR(theme,ntotal,ntheta,contrastLevel,loadData,format,dir)
    if nargin < 4
        dir = '.';
        if nargin < 3
            format ='fig';
            if nargin < 2
                loadData = false;
            end
        end
    end
    dir = [dir,'/'];
    if ~isempty(format)
        if strcmp(format,'psc2')
            printDriver = ['-de',format];
            format = 'eps';
        else
            printDriver = ['-d',format];
        end
        dpi = '-r300';
    end
    if ~loadData
        for i=1:contrastLevel
            disp(['loading c',num2str(i)]);
            eps = 125*2^(i-1);
            conLabel(i) = {[num2str(eps/10,'%.1f'),'%']};
            DIR = [theme,'/',num2str(eps,'%04d')];
            [nP(i),tC(i),aC] = readVsFR(DIR,ntheta,ntotal,ndperiod);
        end
        save([theme,'VsFR.mat'],'tC','nP','aC');
    else 
        load([theme,'VsFR.mat'],'tC','nP','aC');
    end

end
