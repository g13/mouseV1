function neuronlist = spatialFreqTC(fdr,sfLabel,lgnfile,neuronlist,ranking,format,contrast,outputfdr)
	pPosition = [18, 180, 800,700];
	if nargin <8
		outputfdr = fdr{1};
		if nargin < 7
			contrast = 4;
    	    if nargin < 6
				format = '';
                if nargin < 5
                    ranking = {};
                    if nargin < 4
    	                neuronlist = [];
                    end
                end
    	    end
		end
	end
    popOnly = false;
    thres = 0.3;
	kE = 10;
	kI = 2;
	nSF = length(fdr);
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
set(0,'DefaultAxesFontSize',FontSize)
	
    load(lgnfile);
    theta = [reshape(etheta,[p.nv1e,1]);reshape(itheta,[p.nv1i,1])];
    assert(sum(theta>=pi) == 0);
%     theta(theta>=pi) = theta(theta>=pi)-pi;
%     theta(theta<0) = theta(theta<0)+pi;

    ieps = contrast;
    eps = num2str(125*2.^(ieps'-1),'%04d');
    for j=contrast
		d(j) = readData(fdr,eps(j,:),nSF,p.nv1,theta);
        if j==1
            d=repmat(d,[contrast(end),1]);
        end
    % data nSF x p.nv1 x (1,2,3)
          % 1 - prescribed prefered angle
          % 2 - prefered angle by the largest response of simulation result over different SFs
          % 3 - prefered angles by results of different SFs
    end
if ~popOnly 
    if isempty(neuronlist)
        rng(12345);
        neuronlist(1:kE) = randi(p.nv1e,kE,1);
        neuronlist((1:kI)+kE) = p.nv1e + randi(p.nv1i,kI,1);
    end
    LineScheme = {':','-.','--','-',':'};
    % individual e.g.
    titles = {'prescribed angle','max prefered','prefered as is','max input','input as is'};
    for k=1:length(neuronlist)
        i = neuronlist(k);
        hIndividual = figure;
        for itype = 1:5
            subplot(2,3,itype);
            hold on
            for jk = 1:length(contrast)
                j = contrast(jk);
                plot(1:nSF,d(j).fr(:,i,itype),LineScheme{j},'Color','b');
            end
            set(gca,'XTick',1:nSF);
            set(gca,'XTickLabel',sfLabel);
			ylabel('Response at picked Orientation');
			xlabel('SF (10^{-2}cpd)');
			title(titles(itype));
        end
        subplot(2,3,6);
        plot(1:nSF,d(j).cv(:,i),LineScheme{j},'Color','m');
        set(gca,'XTick',1:nSF);
        set(gca,'XTickLabel',sfLabel);
        ylabel('CV');
        xlabel('SF (10^{-2}cpd)');
        if i > p.nv1e
            title({[p.Itypes{p.typeI(i-p.nv1e)}],'CV at diff. SF'});
        else
            title({p.Etypes{p.typeE(i)},'CV at diff. SF'});
        end
        if ~isempty(format)
			set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    	    print(hIndividual,[outputfdr,'/',num2str(i),'_sf_eg.',format],printDriver,dpi);
        end
    end
end
    % population distribution 
    hPop = figure;
    for j = contrast
        subplot(contrast(end),2,2*(j-1) + 1)
        hist(d(j).prefSF,1:nSF);
        set(gca,'XTick',1:nSF);
        set(gca,'XTickLabel',sfLabel);
		xlabel('SF (10^{-2}cpd)');
		title(['pref. SF at ', num2str(j),'th contrast']);

        subplot(contrast(end),2,2*(j-1) + 2)
		increment = (0:nSF:((p.nv1-1)*nSF))';
		ind = d(j).prefSF + increment;
		pick = d(j).pkrate(ind) > d(j).br(ind) + thres;

        hist(d(j).prefSF(pick),1:nSF);
        set(gca,'XTick',1:nSF);
        set(gca,'XTickLabel',sfLabel);
		xlabel('SF (10^{-2}cpd)');
		title(['active only, ',num2str(100*sum(pick)/p.nv1,'%4.4g'), '%']);
    end
    if ~isempty(format)
		set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hPop,[outputfdr,'/sf_pop.',format],printDriver,dpi);
    end

	hType = figure;
	ntE = p.ntypeE;
	ntI = p.ntypeI;
	ylabelsE = strcat(p.Etypes,' Exc');
	nsubE = p.nSubTypeE;
	ylabelsI = strcat(p.Itypes,' Inh');
	nsubI = p.nSubTypeI;
    for j = contrast
        increment = (0:nSF:((p.nv1-1)*nSF))';
        ind = d(j).prefSF + increment;
        pick0 = (d(j).pkrate(ind) > d(j).br(ind) + thres);
		for i = 1:ntE
            subplot(contrast(end),ntE + ntI,(ntE+ntI)*(j-1) + i)
			pick = pick0 & [p.eSubregion == nsubE(i); false(p.nv1i,1)];
			if nsubE(i) == 2 && i == 2
				pick(p.ORF) = false;
			end
			if nsubE(i) == 2 && i == 3
				pick(~p.ORF) = false;
			end

        	hist(d(j).prefSF(pick),1:nSF);
        	set(gca,'XTick',1:nSF);
        	set(gca,'XTickLabel',sfLabel);
			xlabel('SF (10^{-2}cpd)');
			ylabel(ylabelsE{i});
            if i == 1
                title(['prefered SF at ', num2str(j),'th contrast (active neurons only)']);
            end
		end

		for i = 1:ntI
			subplot(contrast(end),ntE + ntI,(ntE+ntI)*(j-1) + ntE + i)
            pick = pick0 & [false(p.nv1e,1); p.iSubregion == nsubI(i)];  
        	hist(d(j).prefSF(pick),1:nSF);
        	set(gca,'XTick',1:nSF);
        	set(gca,'XTickLabel',sfLabel);
			xlabel('SF (10^{-2}cpd)');
			ylabel(ylabelsI{i});
		end
    end
    if ~isempty(format)
% 		set(hType, 'PaperSize',[900 1200],'PaperUnits','centimeters','PaperPosition',[1,1,899,1199])
%         set(hType, 'Units','normalized','OuterPosition',[1,1,1199,899],'Position',[0 0 1.0 1.0]);
%         set(hType, 'PaperType','b3');
		set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
        print(hType,[outputfdr,'/sf_types.',format],printDriver,dpi);
    end


    hOriShift = figure;
    for j = contrast
        for i = 1:nSF
            subplot(contrast(end),nSF,(j-1)*nSF+i)
            hist((d(j).prA(i,:)-d(j).prA(mod(i,nSF)+1,:))*180/pi,10);
            title(['SF',num2str(i),'-SF',num2str(mod(i,nSF)+1),' at ',num2str(j),'th contrast']);
			xlabel('orientation shift (degree)');
        end
    end
    
    if ~isempty(format)
		set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
        print(hOriShift,[outputfdr,'/sf_OriShift.',format],printDriver,dpi);
    end
end
% read in prefered orientation data
function d = readData(fdr,eps,nSF,n,theta)

	% read adg_INPUT for SF, ntheta
    filename = [fdr{1},'/','adg_INPUT'];
    fid = fopen(filename);
    for i=1:5
        fgetl(fid);
    end
    % line 6 
    line = regexprep(fgetl(fid),'\s',' ');
    data = textscan(line,'%f');
    data = data{1};
    ndperiod = data(7);
    tfinal = data(2);

    for i=1:14
        fgetl(fid);
    end
    % line 21
    line = regexprep(fgetl(fid),'\s',' ');
    data = textscan(line,'%f');
    data = data{1};
    %data = data(2:end);
    period = 1/data(1);
    ntheta = data(3);

    for i=1:5
        fgetl(fid);
    end
    % line 27
    line = regexprep(fgetl(fid),'\s',' ');
    data = textscan(line,'%f');
    data = data{1};
    rstart = data(4);
	fclose(fid);
    disp(['period = ',num2str(period),' sampled ',num2str(ndperiod),' times, ntheta = ',num2str(ntheta),'record length = ',num2str(tfinal-rstart),'s']);
    d = struct();
    
	% read cv.dat for extra-cellular properties
    itheta = round(theta*ntheta/pi+1);  % prescribed angles (drifting-grating direction, not bar orientation)
    cv = zeros(n,nSF);
    sc = cv;
    pkrate =cv;
    prA = cv;
    br = cv;
    infrate = cv;
    for i=1:nSF
        DIR = [fdr{i},'/',eps,'/cv.dat'];
        fid = fopen(DIR);
        data = fread(fid, [n, 16],'double');
        cv(:,i) = data(:,1);
        sc(:,i) = data(:,2);
        pkrate(:,i) = data(:,5);
        prA(:,i) = data(:,7);  % actual prefered angle in simulation
%         ei(:,i) = data(:,8);
        br(:,i) = data(:,9);
        priA(:,i) = data(:,12);
        ipriA0 = data(:,15);
        fclose(fid);

        DIR = [fdr{i},'/',eps,'/cycles.dat'];
        fid = fopen([fdr{i},'/',eps, '/frate.dat']);
        data = fread(fid, [n, 2*ntheta],'double');
        fclose(fid);
        for j = 1:n
            infrate(j,i) = data(j,ipriA0(j));
        end
    end
    d.cv = cv';
    d.sc = sc';
    d.pkrate = pkrate';
    d.infrate = infrate';
    d.prA = prA';
    d.br = br';
    d.priA = priA';
	[~, ind] = max(pkrate,[],2);
	d.prefSF = ind;
    increment = (0:nSF:(nSF*(n-1)))';
	ind = ind + increment;
    iprA = round(d.prA*ntheta/pi+1); 
	imaxPrA = iprA(ind);

	[~, ind] = max(infrate,[],2);
	d.prefinSF = ind;
	ind = ind + increment;
    ipriA = round(d.priA*ntheta/pi+1); 
	imaxPriA = ipriA(ind);
    
    % read data of different angles
    increment = (0:2*ntheta:(2*ntheta*(n-1)))';
    temp_imaxPrA = imaxPrA + increment;
    temp_itheta = itheta + increment;
    temp_imaxPriA = imaxPriA + increment;
%     stifr = frate
    d.fr = zeros(nSF,n,5);
    for i=1:nSF
        fid = fopen([fdr{i},'/',eps, '/frate.dat']);
        frate = fread(fid, [n, 2*ntheta],'double');
		frate = reshape(frate',[n*2*ntheta,1]); % transform theta to rows and reshape to vector
%         stifr = fread(fid, [n, 2*ntheta],'double');
        fclose(fid);
        d.fr(i,:,1) = frate(temp_imaxPrA)';
        d.fr(i,:,2) = frate(temp_itheta)';
        d.fr(i,:,4) = frate(temp_imaxPriA)';
    end
    d.fr(:,:,3) = d.pkrate;
    d.fr(:,:,5) = d.infrate;
end
