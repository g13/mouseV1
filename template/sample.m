function sample(theme,lgnfile,processed,all)
    twindow = 1;
    t0 = 0;
    dt = 0.001;
    rtstart = 6-twindow;
    t = rtstart + (t0:dt:twindow);
    tshift = 0;
    %spont = false;
    spont = true;
    eps = 1;
    ntotal = 10800;
    ntheta = 12;
    load(lgnfile);
    nd = [p.enormDistance; p.inormDistance];
    freq = 4;
    period = 1/freq;
    nprint = 20;
    format = 'fig';

    filepath = [theme,'/',num2str(eps*1000,'%0.4d'), '/', 'cv.dat'];
    fid = fopen(filepath);
    Data = fread(fid, [ntotal, 10],'double');
    prA = Data(:,7);
    dtheta = pi/ntheta;
    prA = round(prA/dtheta);
    fclose(fid);
    if all
        [~, sampleID] = plotindividual(theme,lgnfile,[],'',4,12,true,theme,{},0.5,false,true,true);    
        nd = nd(sampleID);
        disp('neuronlist acquired');
        nsample = length(sampleID);
        for i=1:nsample
            disp(['neuron #',num2str(sampleID(i),'%i')]);
            angle = prA(sampleID(i));
            nlgn = nLGN(sampleID(i),1);
            plotsample(theme,eps,angle,sampleID(i),sampleID(i),ntotal,t,nlgn,period,format,processed,lgnfile,tshift,nd(i));
        end
    else
        fid = fopen([theme,'/sample_id.dat']);
        sampleID = fread(fid,[1,inf], 'int');
        nsample = length(sampleID);
        disp([num2str(nsample)]);
        fclose(fid);
        nd = nd(sampleID);
        for i=1:nsample
            disp(['neuron #',num2str(sampleID(i),'%i')]);
            angle = prA(sampleID(i));
            nlgn = nLGN(sampleID(i),1);
            plotsample(theme,eps,angle,i,sampleID(i),nsample,t,nlgn,period,format,processed,lgnfile,tshift,nd(i));
        end
    end
%     

    %if spont  && eps ==0.125 
    %    angle = ntheta; % spont.
    %    nrange = round(linspace(1,nsample,nprint));
    %else
    %    angle = 3;
    %    nrange = find(prA(ID)== angle);
    %    pick = round(linspace(1,length(nrange),nprint));
    %    nrange = nrange(pick);
    %end
%%        nrange = round(linspace(1,nsample,nprint));
%%        angle = 3;
end
