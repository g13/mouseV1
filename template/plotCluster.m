function [] = plotCluster(theme,lgn,ClusterFile,dataRead,format,contrastLevel,ntheta,outputfdr,thres)
% position =[500,150,1000,650];set(0, 'OuterPosition', position);
%position = [50,50,1400,900];set(0, 'DefaultFigurePosition', position);
    pPosition = [0, 0, 1280, 720];
    if nargin < 9
        thres = 0.3; % gauge for active
        if nargin < 8
            outputfdr = theme;
            if nargin < 7
                ntheta = 12;      
                if nargin < 6
                    contrastLevel = 4;
                    if nargin < 5                  
                        format = '';
                    end
                end
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
        dpi = '-r100';
    end
FontSize = 10;
legendFontSizeOffset = 2;
set(0,'DefaultAxesFontSize',FontSize)
% dontUseFit = true;
dontUseFit = false;
tcReady = false;
meanTimeLine = true;
% meanTimeLine = false;
current = true;
% current = false;
rt = 20;
ndperiod = 25;

loadSubstract = false;
polarplot = 0;
pOSI = true;
rng('shuffle');

load(lgn);
load(ClusterFile);
nLGN = nLGN(:,1);
nCluster = p.nv1;
%     figure;
%     subplot(2,2,1);
%     hist(postE); title('postE');
%     subplot(2,2,2);
%     hist(postI); title('postI');
%     subplot(2,2,3);
%     hist(preE); title('preE');
%     subplot(2,2,4);
%     hist(preI); title('preE');
conLabel = cell(contrastLevel,1);
if ~dataRead
    tic;
    for i=1:contrastLevel
        disp(['loading c',num2str(i)]);
        eps = 125*2^(i-1);
        conLabel(i) = {[num2str(eps/10,'%.1f'),'%']};
    %     eps = 250*i;
        DIR = [theme,'/',num2str(eps,'%04d')];
        if i==1
            [tnP,ttC,taC] = readDataAll(DIR,ntheta,nCluster,ndperiod);
            nP = tnP; tC = ttC; aC = taC;
            clear tnP ttC taC
            nP = repmat(nP,contrastLevel,1);
            tC = repmat(tC,contrastLevel,1);
            aC = repmat(oC,contrastLevel,1);
        else
            [nP(i),tC(i),aC(i)] = readDataAll(DIR,ntheta,nCluster,ndperiod);
        end
    end
    lgnmax_e = max(nLGN(nP(1).ei>0.5));
    lgnmin_e = min(nLGN(nP(1).ei>0.5));
    lgnmax_i = max(nLGN(nP(1).ei<0.5));
    lgnmin_i = min(nLGN(nP(1).ei<0.5));
    lowestLGN = min(lgnmin_e,lgnmin_i);
    highestLGN = max(lgnmax_e,lgnmax_i);
    [~,preE,preI] = connections(theme,lgn,nCluster,nP(1).ei');
    
    nv1 = size(nP(i).ei,1);
    orthf = @(rmax,sigfsq2,lifted) rmax.*exp(-(90./sigfsq2).^2)+lifted;
    LineScheme = {':','-.','--','-',':'};
    dtheta = 180/ntheta;
    %% Cluster data gather
    CtC.frate = zeros(nCluster,ntheta);
    CtC.gE = zeros(nCluster,ntheta);
    CtC.gI = zeros(nCluster,ntheta);
    CtC.gLGN = zeros(nCluster,ntheta);
    CtC.gEstd = zeros(nCluster,ntheta);
    CtC.gIstd = zeros(nCluster,ntheta);
    CtC.gLGNsc = zeros(nCluster,ntheta);
    CtC = repmat(CtC,contrastLevel,1);

    CpC.frate = zeros(nCluster,ndperiod);
    CpC.gE = zeros(nCluster,ndperiod);
    CpC.gI = zeros(nCluster,ndperiod);
    CpC.gLGN = zeros(nCluster,ndperiod);
    CpC.gEstd = zeros(nCluster,ndperiod);
    CpC.gIstd = zeros(nCluster,ndperiod);
    CpC.gLGNstd = zeros(nCluster,ndperiod);
    CpC = repmat(CpC,contrastLevel,1);

    CipC.frate = zeros(nCluster,ndperiod);
    CipC.gE = zeros(nCluster,ndperiod);
    CipC.gI = zeros(nCluster,ndperiod);
    CipC.gLGN = zeros(nCluster,ndperiod);
    CipC.gEstd = zeros(nCluster,ndperiod);
    CipC.gIstd = zeros(nCluster,ndperiod);
    CipC.gLGNstd = zeros(nCluster,ndperiod);
    CipC = repmat(CipC,contrastLevel,1);

    CnP.iPrA = zeros(nCluster,1);
    CnP.iPriA = zeros(nCluster,1);
    CnP.CV = zeros(nCluster,1);
    CnP = repmat(CnP,contrastLevel,1);

    cCVsortedID = zeros(nCluster,contrastLevel);
    cPKsortedID = zeros(nCluster,contrastLevel);
    
    rtheta = linspace(0,2*pi,2*ntheta+1);
    phase = linspace(0,2*pi,ndperiod+1);
    iphase = 1i*ones(nCluster,1)*phase(1:ndperiod);
    irtheta = 1i*ones(nCluster,1)*rtheta(1:2*ntheta);
    for i=1:nCluster
        Cbr = mean(nP(contrastLevel).br(cluster{i}));
        CnLGN = nearest(mean(nLGN(cluster{i},1)));
    end
    for j=1:contrastLevel
        for i=1:nCluster
            %% circular mean prefered angles
            Angles = (nP(i).indpo(cluster{i})-1)*dtheta/180*pi;
            CnP(j).iPrA(i) = nearest(angle(mean(exp(i*Angles)))/(dtheta/180*pi));
            Angles = (nP(i).indpi(cluster{i})-1)*dtheta/180*pi;
            CnP(j).iPriA(i) = nearest(angle(mean(exp(i*Angles)))/(dtheta/180*pi));

            CtC(j).frate(i,:) = mean(tC(j).frate(cluster{i},:),1);
            CtC(j).gE(i,:) = mean(tC(j).gE(cluster{i},:),1);
            CtC(j).gI(i,:) = mean(tC(j).gI(cluster{i},:),1);
            CtC(j).gLGN(i,:) = mean(tC(j).gLGN(cluster{i},:),1);
            CtC(j).gEstd(i,:) = mean(tC(j).gEstd(cluster{i},:),1);
            CtC(j).gIstd(i,:) = mean(tC(j).gIstd(cluster{i},:),1);
            CtC(j).gLGNsc(i,:) = mean(tC(j).gLGNsc(cluster{i},:),1);

            CpC(j).frate(i,:)=mean(aC(j).spikes(cluster{i},:,CnP(j).iPrA(i))*ndperiod/rt,1);
            CpC(j).gE(i,:) = mean(aC(j).gE(cluster{i},:,CnP(j).iPrA(i)),1);
            CpC(j).gI(i,:) = mean(aC(j).gI(cluster{i},:,CnP(j).iPrA(i)),1);
            CpC(j).gLGN(i,:) = mean(aC(j).gLGN(cluster{i},:,CnP(j).iPrA(i)),1);
            CpC(j).gEstd(i,:) = mean(aC(j).gEstd(cluster{i},:,CnP(j).iPrA(i)),1);
            CpC(j).gIstd(i,:) = mean(aC(j).gIstd(cluster{i},:,CnP(j).iPrA(i)),1);
            CpC(j).gLGNstd(i,:) = mean(aC(j).gLGNsc(cluster{i},:,CnP(j).iPrA(i)),1);

            CipC(j).frate(i,:)=mean(aC(j).spikes(cluster{i},:,CnP(j).iPrA(i))*ndperiod/rt,1);
            CipC(j).gE(i,:) = mean(aC(j).gE(cluster{i},:,CnP(j).iPriA(i)),1);
            CipC(j).gI(i,:) = mean(aC(j).gI(cluster{i},:,CnP(j).iPriA(i)),1);
            CipC(j).gLGN(i,:) = mean(aC(j).gLGN(cluster{i},:,CnP(j).iPriA(i)),1);
            CipC(j).gEstd(i,:) = mean(aC(j).gEstd(cluster{i},:,CnP(j).iPriA(i)),1);
            CipC(j).gIstd(i,:) = mean(aC(j).gIstd(cluster{i},:,CnP(j).iPriA(i)),1);
            CipC(j).gLGNstd(i,:) = mean(aC(j).gLGNsc(cluster{i},:,CnP(j).iPriA(i)),1);
        end 
        %% CV & sortIDs
        CnP(j).PK = max(CtC(j).frate,[],2);
        CnP(j).MR = mean(CtC(j).frate,2);
        CnP(j).CV = 1-abs(sum(exp(2*irtheta).*[CtC(j).frate,CtC(j).frate],2)./(ntheta*CnP(j).MR));
        CnP(j).SC = 2*abs(sum(exp(iphase).*CipC(j).frate,2)./sum(CipC(j).frate,2));
        [~,cCVsortedID] = sort(CnP(j).CV);
        [~,cPKsortedID] = sort(CnP(j).PK);
    end
    [~,cSizeSortedID] = sort(clusterSize);
    toc;
    save([theme,'-tcData-Cx',num2str(contrastLevel),'.mat'],'CtC','CpC','CipC','CnP','cCVsortedID','cPKsortedID','cSizeSortedID','Cbr','CnLGN');
else
    load([theme,'-tcData-Cx',num2str(contrastLevel),'.mat']);
end
for i=1:contrastLevel
    if tcReady
        rmax = zeros(contrastLevel,p.nv1);
        smax = rmax;
        sigfsq2 = rmax;
        lifted = rmax;
        if loadSubstract
            if exist([theme,'/',theme,'-',num2str(i),'fitted_substracted.mat'],'file');
                fitted = load([theme,'/',theme,'-',num2str(i),'fitted_substracted.mat']);
                rmax(i,:) = fitted.rmax;
                smax(i,:) = fitted.smax;
                sigfsq2(i,:) = fitted.sigfsq2;
                lifted(i,:) = fitted.lifted;
                tcReady = true;
            else
                disp(['contrast ',num2str(i),' fitted curve not ready']);
                tcReady = false;
            end
        else
            if exist([theme,'/',theme,'-',num2str(i),'fitted.mat'],'file');
                fitted = load([theme,'/',theme,'-',num2str(i),'fitted.mat']);
                rmax(i,:) = fitted.rmax;
                smax(i,:) = fitted.smax;
                sigfsq2(i,:) = fitted.sigfsq2;
                lifted(i,:) = fitted.lifted;
                tcReady = true;
            else
                disp(['contrast ',num2str(i),' fitted curve not ready']);
                tcReady = false;
            end
        end
    end
end
disp('data loaded');
if ~statsOnly
    level = contrastLevel;
    fired = CnP(level).PK(cCVsortedID(:,level)) > max(Cbr(cCVsortedID(:,level)),thres);
    sn = 5;
    clusterList = zeros(sn*3+3, 1);
    ranking = cell(length(clusterList),1);

    otherType = ' Simple';
    j = 0;
    type0 = CnP(level).SC > 1.0; 
    types = type0(cCVsortedID(:,level));
    candidates = find(fired & types);
    ncand = length(candidates);
    if ncand>=sn-2
        clusterList(j + (1:sn-2)) = cCVsortedID(candidates([1, floor(ncand/2), ncand]),level);
        ranking(j + (1:sn-2)) = {'smallest CV','median CV','largest CV'};
    end
    types = type0(cPKsortedID(:,level));
    candidates = find(types);
    ncand = length(candidates);
    if ncand>=2
        clusterList(j + (sn-1:sn)) = cPKsortedID(candidates([floor(ncand/2),ncand]),level);
        ranking(j + (sn-1:sn)) = {'median PK','largest PK'};
    end
        
    otherType = ' Complex';
    j = sn*1;
    type0 = CnP(level).SC < 1.0; 
    types = type0(cCVsortedID(:,level));
    candidates = find(fired & types);
    ncand = length(candidates);
    if ncand>=sn-2
        clusterList(j + (1:sn-2)) = cCVsortedID(candidates([1, floor(ncand/2), ncand]),level);
        ranking(j + (1:sn-2)) = strcat({'smallest CV','median CV','largest CV'},otherType);
    end
    disp([num2str(ncand),' active ',otherType,' clusters are found']);
    if ncand >= 2
        types = type0(cPKsortedID(:,level));
        candidates = find(types);    
        ncand = length(candidates);
        if ncand>=2
            clusterList(j + (sn-1:sn)) = cPKsortedID(candidates([floor(ncand/2),ncand]),level);
            ranking(j + (sn-1:sn)) = strcat({'median PK','largest PK'},otherType);
        end
    end

    %%
    clusterSize = clusterSize(1:nCluster);
    clusterList(2*sn+(1:3)) = [cSizeSortedID(1),cSizeSortedID(round(nCluster/2)),cSizeSortedID(nCluster)];
    ranking(2*sn+(1:3)) = {'smallest Size','median Size','largest Size'}; 
    pickedNeuron = clusterList ~= 0;
    clusterList =  clusterList(pickedNeuron);
    ranking = ranking(pickedNeuron);
if sum(clusterList) == 0
    disp('No single neuron match the filter');
else
    tic;
    FWHM = 2*sqrt(2*log(2));
    HWHM = sqrt(2*log(2));
    bound = p.bound;
    rfx = @(x0,a,b,theta,t,ra,rb) x0 + (ra.*a.*cos(t).*cos(theta)-rb.*b.*sin(t).*sin(theta));
    rfy = @(y0,a,b,theta,t,ra,rb) y0 + (ra.*a.*cos(t).*sin(theta)+rb.*b.*sin(t).*cos(theta));
    ttt = 0:0.1:2*pi;
    for nn = 1:length(clusterList)
        k = clusterList(nn);
        
        ipA = CnP(contrastLevel).ipriA(k);
        inputAngle = (ipA-1)*dtheta;

        if ipA < ntheta/2+1
            half = [(ipA+ntheta/2):ntheta,1:ipA+ntheta/2];
        else 
            half = [(ipA-ntheta/2):ntheta,1:(ipA-ntheta/2)];
        end
        if inputAngle >= pi/2
            pA = inputAngle - pi/2;
        else
            pA = inputAngle + pi/2;
        end
        grating_ipA = round(pA*180/pi/dtheta)+1;
        thetas = ((grating_ipA-ntheta/2-1):(grating_ipA+ntheta/2-1))*dtheta;
    
        hIndividual = figure;
        subplot(2,3,1); 
        for i=1:contrastLevel
    
            rho = [CtC(i).frate(k,:),CtC(i).frate(k,1)];
            if polarplot
                rtheta = linspace(0,2*pi,2*ntheta+1);
                polar(rtheta,rho,LineScheme{i});
            else
                plot(thetas,rho(half),'Color','k','LineStyle',LineScheme{i});
            end
            hold on
        end
        if polarplot
            polar(rtheta,Cbr(k)*ones(1,17),':r');
        else
            plot(thetas(end),Cbr(k),'*r');
            xlim([thetas(1),thetas(end)]);
        end
        ylabel('firing rate');
    
        if tcReady && ~dontUseFit
            for il = 1:contrastLevel
                s = (smax(il,k)-90):1:(smax(il,k)+90);
                r = rmax(il,k)*exp(-((s-smax(il,k))./sigfsq2(il,k)).^2)+lifted(il,k);
                plot(s,r,LineScheme{il},'Color','g','LineWidth',2);
            end
        end
        yy = ylim();
        ylim([0,yy(2)*1.2]);
        legend(conLabel,'FontSize',FontSize-legendFontSizeOffset);
    
        cvL = zeros(contrastLevel,1);
        cvLnB = zeros(contrastLevel,1);
        for i=1:contrastLevel
            cvL(i) = CnP(i).CV(k);
            cvLnB(i) = CnP(i).CVNoBack(k);
        end
        cvLevel = strjoin(cellstr(num2str(cvL,'%1.1f'))'); 
        cvNBLevel = strjoin(cellstr(num2str(cvLnB,'%1.1f'))');
    %     for i=1:contrastLevel
    %         cvLevel = [cvLevel,' ', num2str(CnP(i).CV(k),'%1.2f')];
    %         cvNBLevel = [cvNBLevel,' ', num2str(CnP(i).CVNoBack(k),'%1.2f')];
    %     end
        % title({['frate',', CV =', cvLevel],['stimulated CV =',cvNBLevel],['\theta = ', num2str(ptheta(k)*180/pi),'^o']});
        title({['CV =', cvLevel],['evoked CV =',cvNBLevel]});
    
       subplot(2,3,2); hold on;
        for i=1:contrastLevel
            plot(thetas, CtC(i).gLGN(k,half),'LineStyle',LineScheme{i},'Color','g');
            plot(thetas, CtC(i).gLGN(k,half).*(1+CtC(i).gLGNsc(k,half)),'LineStyle',LineScheme{i},'Color','k');
            plot(thetas, CtC(i).gE(k,half),'LineStyle',LineScheme{i},'Color','r');
            plot(thetas, CtC(i).gE(k,half)+ CtC(i).gEstd(k,half),'LineStyle',LineScheme{i},'Color','b');
            plot(thetas, CtC(i).gEn(k,half),'.r');
        end
        xlim([thetas(1),thetas(end)]);
        if CnP(1).ei(k) > 0.5
            neuron_alias = ['Cluster-',ranking{nn}];
        else
            neuron_alias = ['Cluster-',ranking{nn}];
        end
        title({['No.', num2str(k)],neuron_alias} );
        yy = ylim();
        ylim([0,yy(2)*1.3]);
        ylabel('Conductance');
        xlabel(['relative angle to ', conLabel{contrastLevel},' input']);
        legend({'gLGN','gLGN-F1','gE','gEstd','gEn'},'FontSize',FontSize-legendFontSizeOffset);
        
        subplot(2,3,3);hold on;
        for i=1:contrastLevel
            mmm = max(CtC(i).gE(k,half));
            plot(thetas, CtC(i).gE(k,half)./mmm,'LineStyle',LineScheme{i},'Color','r');
            
            mmm = max(CtC(i).gI(k,half));
            plot(thetas, CtC(i).gI(k,half)./mmm,'LineStyle',LineScheme{i},'Color','b');
            
            mmm = max(CtC(i).gLGNsc(k,half));
            plot(thetas, CtC(i).gLGNsc(k,half)./mmm,'LineStyle',LineScheme{i},'Color','k');
    
            mmm = max(CtC(i).gLGN(k,half));
            plot(thetas, CtC(i).gLGN(k,half)./mmm,'LineStyle',LineScheme{i},'Color','g');
        end
        title({['\lambda ~',num2str((nLGN(k)-lgnmin_e)/(lgnmax_e-lgnmin_e),'%1.2f'),', nLGN=',num2str(nLGN(k))],...
                ['preE: ',num2str(preE(k)),', preI: ',num2str(preI(k))]});
    %     [',distance to pinwheel = ', num2str(CnP(1).dist(k))];
        ylabel('Normalized g');
        xlim([thetas(1),thetas(end)]);
    
        subplot(2,3,4); hold on;
        for i=1:contrastLevel
            plot(thetas, CtC(i).Veff(k,half),'LineStyle',LineScheme{i},'Color','b');
            plot(thetas, CtC(i).Veff(k,half).*(1+CtC(i).Veffsc(k,half)),'LineStyle',LineScheme{i},'Color','r');
        end
        xlim([thetas(1),thetas(end)]);
        ylabel('V_{eff}');
        yy = ylim();
        ylim([0,yy(2)*1.2]);
        legend({'Veff','Veff-F1'},'FontSize',FontSize-legendFontSizeOffset);
    
        OSI = zeros(contrastLevel,1);
        OSInB = OSI;
        prefinTheta = OSI;
        prefoutTheta = OSI;
        if tcReady && ~dontUseFit
            for i = 1:contrastLevel
                rp = rmax(i,k)+lifted(i,k);
                ro = orthf(rmax(i,k),sigfsq2(i,k),lifted(i,k));
                OSI(i) = (rp-ro)/(ro+rp);
                OSInB(i) = (rp-ro)/(ro+rp- 2*Cbr(k));
            end
            pick = smax(i,:)>=90;
            prefoutTheta(pick) = smax(i,pick)-90;
            prefoutTheta(~pick) = smax(i,~pick)+90;
        else
            for i = 1:contrastLevel
                iprefTheta = CnP(i).iPrA(k);
                prefoutTheta(i) = iprefTheta*dtheta/180*pi;
                iorthTheta = mod(iprefTheta+ntheta/2-1, ntheta)+1;
                OSI(i) = (CtC(i).frate(k,iprefTheta) - CtC(i).frate(k,iorthTheta))/...
                    (CtC(i).frate(k,iprefTheta) + CtC(i).frate(k,iorthTheta));
                OSInB(i) = (CtC(i).frate(k,iprefTheta) - CtC(i).frate(k,iorthTheta))/...
                    (CtC(i).frate(k,iprefTheta) + CtC(i).frate(k,iorthTheta) - 2*Cbr(k));
                %oangles = [oangles,num2str(prefoutTheta,'%3.0f'),'^{o} '];
            end
            pick = prefoutTheta>=90;
            prefoutTheta(pick) = prefoutTheta(pick)-90;
            prefoutTheta(~pick) = prefoutTheta(~pick)+90;
        end
        oangles = strjoin(cellstr(strcat(num2str(prefoutTheta,'%3.0f'),'^{o}'))');
        for i = 1:contrastLevel
            prefinTheta(i) = CnP(i).priA(k)*180/pi;
        end
        pick = prefinTheta>=90;
        prefinTheta(pick) = prefinTheta(pick)-90;
        prefinTheta(~pick) = prefinTheta(~pick)+90;
    
        iangles = strjoin(cellstr(strcat(num2str(prefinTheta,'%3.0f'),'^{o}'))');
    
        OSI(isnan(OSI)) = 0;
        OSInB(isnan(OSInB)) = 0;
        sOSI =  strjoin(cellstr(num2str(OSI,'%1.1f'))');
        sOSInB =  strjoin(cellstr(num2str(OSInB,'%1.1f'))');
        title({['input angle :', iangles],['output angle :', oangles]});
    
        subplot(2,3,5); hold on;
    	subregion = length(nSubLGN{k});
        bounda = bound*ones(1,subregion);
        subpick = 1:subregion;
        if k <= p.nv1e
            theta = etheta(k);
    		sigma = p.esigma(k,subpick);
        	sigmb = sigma.*p.eAspectRatio(k,subpick);
            boundb = bounda.*p.eAspectRatio(k,subpick);
        else
            theta = itheta(k-p.nv1e);
    		sigma = p.isigma(k-p.nv1e,subpick);
        	sigmb = sigma.*p.iAspectRatio(k-p.nv1e,subpick);
            boundb = bounda.*p.iAspectRatio(k-p.nv1e,subpick);
        end
        if theta >=pi/2
            gtheta = theta - pi/2;
        else
            gtheta = theta + pi/2;
        end
    
        peak = zeros(subregion,2);
        for j = 1:subregion
            if nSubLGN{k}(j) > 0
                peak(j,:) = mean(LGNpos(v1Map{k,j},:),1)*180/pi;
            else
                peak(j,:) = mean(LGNpos(-v1Map{k,j},:),1)*180/pi;
            end
        end
        ymax = -inf;
        ymin = inf;
        xmax = ymax;
        xmin = ymin;
        
        for j = 1:subregion
            xxx0 = rfx(peak(j,1)/180*pi,sigmb(j),sigma(j),gtheta,ttt,boundb(j),bounda(j))*180/pi;         
            yyy0 = rfy(peak(j,2)/180*pi,sigmb(j),sigma(j),gtheta,ttt,boundb(j),bounda(j))*180/pi;
            xxx1 = rfx(peak(j,1)/180*pi,sigmb(j),sigma(j),gtheta,ttt,HWHM,HWHM)*180/pi;         
            yyy1 = rfy(peak(j,2)/180*pi,sigmb(j),sigma(j),gtheta,ttt,HWHM,HWHM)*180/pi;
            xxx2 = rfx(peak(j,1)/180*pi,sigmb(j),sigma(j),gtheta,ttt,FWHM,FWHM)*180/pi;         
            yyy2 = rfy(peak(j,2)/180*pi,sigmb(j),sigma(j),gtheta,ttt,FWHM,FWHM)*180/pi;
            ymin = min(min(yyy2),ymin);
            ymax = max(max(yyy2),ymax);
            xmin = min(min(xxx2),xmin);
            xmax = max(max(xxx2),xmax);
            x = [xmin,xmax];
            if nSubLGN{k}(j) > 0   
                plot(xxx2,yyy2,'--r','LineWidth',1);
                plot(xxx1,yyy1,':r','LineWidth',1); 
                plot(xxx0,yyy0,'-r','LineWidth',1); 
                plot(peak(j,1),peak(j,2),'sr');
                %plot3(LGNpos(v1Map{k,j},1)*180/pi, LGNpos(v1Map{k,j},2)*180/pi,lgnStrength{k,j}*5,'*r','MarkerSize',4.0);
                plot(LGNpos(v1Map{k,j},1)*180/pi, LGNpos(v1Map{k,j},2)*180/pi,'*r','MarkerSize',4.0);
                plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':r');
            else
                plot(xxx2,yyy2,'--b','LineWidth',1);
                plot(xxx1,yyy1,':b','LineWidth',1); 
                plot(xxx0,yyy0,'-b','LineWidth',1); 
                plot(peak(j,1),peak(j,2),'sb');
                %plot3(LGNpos(-v1Map{k,j},1)*180/pi, LGNpos(-v1Map{k,j},2)*180/pi,lgnStrength{k,j}*5,'ob','MarkerSize',4.0);
                plot(LGNpos(-v1Map{k,j},1)*180/pi, LGNpos(-v1Map{k,j},2)*180/pi,'ob','MarkerSize',4.0);
                plot(x,tan(gtheta)*(x-peak(j,1))+peak(j,2),':b');
            end
        end
        grid on
        axis equal
        ylim([ymin,ymax]);
        xlim([xmin,xmax]);
        dxy = p.rlgn*180/pi;
        edge = 1;
        x = (xmin-edge*dxy):(dxy/20):(xmax+edge*dxy);
        y = (ymin-edge*dxy):(dxy/20):(ymax+edge*dxy);   
    
    %     zmax = 0;
        for j = 1:subregion
           if nSubLGN{k}(j) > 0
                [X, Y, Z] = multiLGNspatialKernel(LGNpos(v1Map{k,j},1)*180/pi, LGNpos(v1Map{k,j},2)*180/pi,x,y,nSubLGN{k}(j),lgnStrength{k,j},200);
                contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','m','LineWidth',2);
                newZmax = max(max(Z));
                contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','m','LineWidth',2);
           else
                [X, Y, Z] = multiLGNspatialKernel(LGNpos(-v1Map{k,j},1)*180/pi, LGNpos(-v1Map{k,j},2)*180/pi,x,y,-nSubLGN{k}(j),lgnStrength{k,j},200);
                contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','g','LineWidth',2);
                newZmax = max(max(Z));
                contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','g','LineWidth',2);
           end
    %        mesh(X,Y,Z,'FaceColor','none');
    %        if  newZmax > zmax
    %            zmax = newZmax;
    %        end
        end
        ylim([y(1),y(end)]);
        xlim([x(1),x(end)]);
        xlabel('degree');
        ylabel('degree');
        %campos([xrange(1),yrange(1),zmax*1.5]);
    
        scL = zeros(contrastLevel,1);
        for i = 1:contrastLevel
            scL(i) = CnP(i).SC(k);
        end
        scLevel = strjoin(cellstr(num2str(scL,'%1.1f'))');
    
        title({['preset grating angle: ',num2str(gtheta*180/pi,'%3.0f'),'^{o}, #sub:',num2str(reshape(nSubLGN{k},[1,subregion]),'%2i')],...
    			['F1/F0 =', scLevel]});
    
        subplot(2,3,6); hold on;
        for i=1:contrastLevel
            errorbar(thetas+5*(i-1), CtC(i).gI(k,half),CtC(i).gIstd(k,half),'LineStyle',LineScheme{i},'Color','b');
            plot(thetas, CtC(i).gIn(k,half),'.b');
        end
        title({['OSI: ',sOSI],['OSInB: ',sOSInB]});
        xlim([thetas(1),thetas(end)]);
        yy = ylim();
        ylim([0,yy(2)*1.2]);
        legend({'gI','gIn'},'FontSize',FontSize-legendFontSizeOffset);
        if ~isempty(format)
            set(gcf,'Renderer','Painters')
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition*1.5);
            print(hIndividual,[outputfdr,'/',num2str(k),'-',neuron_alias,'-tc-',theme,'.',format],printDriver,dpi);
        end
        %% current
        if current
            hCurrent = figure;
            for i=1:contrastLevel
                subplot(2,2,1)
                hold on
                errorbar(thetas+5*(i-1), CtC(i).currEt(k,half),CtC(i).currEtstd(k,half),'LineStyle',LineScheme{i},'Color','g');
                if i == 1, title('LGN current');xlim([thetas(1),thetas(end)]);end
                
                subplot(2,2,2)
                hold on
                errorbar(thetas+5*(i-1), CtC(i).currEn(k,half),CtC(i).currEnstd(k,half),'LineStyle',LineScheme{i},'Color','r');
                errorbar(thetas+5*(i-1), CtC(i).currIn(k,half),CtC(i).currInstd(k,half),'LineStyle',LineScheme{i},'Color','b');
                if i == 1, title('E and I noise current');xlim([thetas(1),thetas(end)]);end
                
                subplot(2,2,3)
                hold on
                errorbar(thetas+5*(i-1), CtC(i).currEc(k,half),CtC(i).currEcstd(k,half),'LineStyle',LineScheme{i},'Color','r');
                if i == 1, title('Exc Cortical current');xlim([thetas(1),thetas(end)]);end
                
                subplot(2,2,4)
                hold on
                errorbar(thetas+5*(i-1), CtC(i).currIc(k,half),CtC(i).currIcstd(k,half),'LineStyle',LineScheme{i},'Color','b');
                if i == 1, title('Inh Cortical current');xlim([thetas(1),thetas(end)]);end
            end
            if ~isempty(format)
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                print(hCurrent,[outputfdr,'/',num2str(k),'-',neuron_alias,'-current-',theme,'.',format],printDriver,dpi);
            end
        end
        %% timeline average
        if meanTimeLine
            hTimeline = figure;
            for i=1:contrastLevel
                
                iprefTheta = CnP(i).iPrA(k);
                subplot(contrastLevel,4,(i-1)*4+1);hold on;
                
                errorbar(1:ndperiod,CpC(i).gLGN(k,:),CpC(i).gLGNstd(k,:),'-g');
                plot(ndperiod+0.6,mean(CpC(i).gLGN(k,:)),'sg');
                errorbar((1:ndperiod)+0.3,oC(i).gLGN(k,:),oC(i).gLGNstd(k,:),':k');
                plot(ndperiod+0.3,mean(oC(i).gLGN(k,:)),'^k');
                if i == 1
                    title(['gLGN, F1/F0: ',num2str(CtC(i).gLGNsc(k,iprefTheta),'%1.2f')]); 
                else
                    title(num2str(CtC(i).gLGNsc(k,iprefTheta),'%1.2f')); 
                end
                ylabel(conLabel{i});
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+2);hold on;
                
                errorbar(1:ndperiod,CpC(i).gI(k,:), CpC(i).gIstd(k,:),'-b');
                plot(ndperiod+0.6,mean(CpC(i).gI(k,:)),'sb');
                errorbar(1:ndperiod,oC(i).gI(k,:), oC(i).gIstd(k,:),':k');
                plot(ndperiod,mean(oC(i).gI(k,:)),'^k');
                if i == 1, title(['gInh, ',neuron_alias]); end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+3); hold on;  
                errorbar(1:ndperiod,CpC(i).gE(k,:), CpC(i).gEstd(k,:),'-r');
                plot(ndperiod+0.6,mean(CpC(i).gE(k,:)),'sr');
                errorbar((1:ndperiod)+0.3,oC(i).gE(k,:), oC(i).gEstd(k,:),':k');
                plot(ndperiod+0.3,mean(oC(i).gE(k,:)),'^k');
                if i == 1
                    title(['gE, F1/F0: ',num2str(CtC(i).gEsc(k,iprefTheta),'%1.2f')]);
                else
                    title(num2str(CtC(i).gEsc(k,iprefTheta),'%1.2f'));
                end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+4); hold on;
                if i==1,title(['FR, ID.', num2str(k)]);end
                plot(1:ndperiod,CpC(i).spikes(k,:)*ndperiod/rt,'-k');
                plot(1:ndperiod,oC(i).spikes(k,:)*ndperiod/rt,':k');
                plot(ndperiod+0.6,mean(CpC(i).spikes(k,:))*ndperiod/rt,'sk');
                plot(ndperiod+0.3,mean(oC(i).spikes(k,:))*ndperiod/rt,'^k');
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                if i==contrastLevel, xlabel('phase');end
            end
            if ~isempty(format)
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                print(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-operiod-',theme,'.',format],printDriver,dpi);
            end
    
            hTimeline = figure;
            for i=1:contrastLevel
                
                iprefTheta = CnP(i).iPriA(k);
                subplot(contrastLevel,4,(i-1)*4+1);hold on;
                
                errorbar(1:ndperiod,CipC(i).gLGN(k,:),CipC(i).gLGNstd(k,:),'-g');
                plot(ndperiod+0.6,mean(CipC(i).gLGN(k,:)),'sg');
                errorbar((1:ndperiod)+0.3,ioC(i).gLGN(k,:),ioC(i).gLGNstd(k,:),':k');
                plot(ndperiod+0.3,mean(ioC(i).gLGN(k,:)),'^k');
                if i == 1
                    title(['gLGN, F1/F0: ',num2str(CtC(i).gLGNsc(k,iprefTheta),'%1.2f')]); 
                else
                    title(num2str(CtC(i).gLGNsc(k,iprefTheta),'%1.2f')); 
                end
                ylabel(conLabel{i});
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+2);hold on;
                
                errorbar(1:ndperiod,CipC(i).gI(k,:), CipC(i).gIstd(k,:),'-b');
                plot(ndperiod+0.6,mean(CipC(i).gI(k,:)),'sb');
                errorbar(1:ndperiod,ioC(i).gI(k,:), ioC(i).gIstd(k,:),':k');
                plot(ndperiod,mean(ioC(i).gI(k,:)),'^k');
                if i == 1, title(['gInh, ',neuron_alias]); end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+3); hold on;  
                errorbar(1:ndperiod,CipC(i).gE(k,:), CipC(i).gEstd(k,:),'-r');
                plot(ndperiod+0.6,mean(CipC(i).gE(k,:)),'sr');
                errorbar((1:ndperiod)+0.3,ioC(i).gE(k,:), ioC(i).gEstd(k,:),':k');
                plot(ndperiod+0.3,mean(ioC(i).gE(k,:)),'^k');
                if i == 1
                    title(['gE, F1/F0: ',num2str(CtC(i).gEsc(k,iprefTheta),'%1.2f')]);
                else
                    title(num2str(CtC(i).gEsc(k,iprefTheta),'%1.2f'));
                end
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                
                subplot(contrastLevel,4,(i-1)*4+4); hold on;
                if i==1,title(['FR, ID.', num2str(k)]);end
                plot(1:ndperiod,CipC(i).spikes(k,:)*ndperiod/rt,'-k');
                plot(1:ndperiod,ioC(i).spikes(k,:)*ndperiod/rt,':k');
                plot(ndperiod+0.6,mean(CipC(i).spikes(k,:))*ndperiod/rt,'sk');
                plot(ndperiod+0.3,mean(ioC(i).spikes(k,:))*ndperiod/rt,'^k');
                ylim([0,inf]);
                xlim([0,ndperiod+2]);
                if i==contrastLevel, xlabel('phase');end
            end
            if ~isempty(format)
                set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
                print(hTimeline,[outputfdr,'/',num2str(k),'-',neuron_alias,'-iperiod-',theme,'.',format],printDriver,dpi);
            end
        end
    end
    disp('individual examples done');
    toc;
end
end

%%
if (~individOnly)
    tic;
    hCVpair = figure;
    hold on; 
    cvLGNmean = zeros(highestLGN-lowestLGN+1,contrastLevel)-1;
    %cvLGNstd = cvLGNmean;
    p1 =contrastLevel-2; p2 = contrastLevel;
    while p1<=0
        p1 = p1 + 1;
    end
    %suptitle(theme);
    for i=lgnmin_e:lgnmax_e
        pick= CnLGN==i & CnP(p1).PK> Cbr & CnP(p2).PK> Cbr & CnP(contrastLevel).PK>thres;
        if ~isempty(pick)
            plot(1-CnP(p1).CV(pick),1-CnP(p2).CV(pick),'o','Color',[0.1+0.9*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
        end
        for j = 1:contrastLevel
            pick= CnLGN==i & CnP(j).PK> Cbr & CnP(contrastLevel).PK>thres;
            if ~isempty(pick)
                cvLGNmean(i-lowestLGN+1,j) = mean(CnP(j).CV(pick));
            end
        end
    end
    plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
    ylabel(['1-CV at ',num2str(12.5*2^(p2-1),'%3.1f'),'% Contrast']);
    xlabel(['1-CV at ',num2str(12.5*2^(p1-1),'%3.1f'),'% Contrast']);
    title('Excitatory');
    
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hCVpair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_CV-Cluster-',theme,'.',format],printDriver,dpi);
    end
    
    hCVnLGN = figure;
    hold on;
    x = (-0.2:(highestLGN-lowestLGN-0.2));
    for i=1:contrastLevel
    %    x = (-0.2:(highestLGN-lowestLGN-0.2))+0.1*(i-1);
        pick = cvLGNmean(:,i)>=0;
        plot(x(pick),1-cvLGNmean(pick,i),'-*','Color',[0.1+0.9*i/contrastLevel,0,0]);
    end
    xlabel('# LGN input');
    ylabel('1-CV');
    title('active Exc 1-CV');
    xlim([lgnmin_e, lgnmax_e]);
    ylim([0,inf]);
    yy = ylim();
    ylim([yy(1),yy(2)*1.2]);
    legend(conLabel);
    
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hCVnLGN,[outputfdr,'/','CVoverLGN-Cluster-',theme,'.',format],printDriver,dpi);
    end
    % subplot(1,2,2); hold on;
    % for i=1:contrastLevel
    %     x = 0:lgnmax+0.2;
    %     errorbar(x,cvLGNmean(:,i,2),cvLGNstd(:,i,2),'*','Color',[0,0,0.1+0.9*i/lgnmax]);
    %     bar(x,cvLGNmean(i,:,2),'FaceColor',[0.1+0.9*i/lgnmax,0,0],'BarWidth',0.1);
    % end
    %% OSI
    hPrefOrth = figure; 
    eoorate = zeros(contrastLevel,nCluster);
    eporate = zeros(contrastLevel,nCluster);
    eoirate = zeros(contrastLevel,nCluster);
    epirate = zeros(contrastLevel,nCluster);
    
    aiprefTheta = zeros(contrastLevel,nCluster);
    aiorthTheta = aiprefTheta;
    for i = 1:contrastLevel
        aiprefTheta(i,:) = CnP(i).iPrA;
        aiorthTheta(i,:) = mod(aiprefTheta(i,:)+ntheta/2-1,ntheta)+1;
    end
    
    
    for j=1:contrastLevel
        nntheta = size(CtC(j).frate,2);
        %pref
        pick = false(nCluster * nntheta,1);
        opPick = aiprefTheta(j,:) + (0:nntheta:((nCluster-1)*nntheta));
        pick(opPick) = true;
        pick = reshape(pick,[nntheta,nCluster])';
        
        eporate(j,:) = CtC(j).frate(pick)';
        
        % orth
        pick = false(nCluster * nntheta,1);
        opPick = aiorthTheta(j,:) + (0:nntheta:((nCluster-1)*nntheta));
        pick(opPick) = true;
        pick = reshape(pick,[nntheta,nCluster])';
        
        eoorate(j,:) = CtC(j).frate(pick)';
    end
    
    
    aiprefTheta_in = zeros(contrastLevel,nCluster);
    aiorthTheta_in = aiprefTheta_in;
    for i = 1:contrastLevel
        aiprefTheta_in(i,:) = CnP(i).iPriA;
        aiorthTheta_in(i,:) = mod(aiprefTheta_in(i,:)+ntheta/2-1,ntheta)+1;
    end
    
    for j=1:contrastLevel
        nntheta = size(CtC(j).frate,2);
        %pref
        pick = false(nCluster * nntheta,1);
        opPick = aiprefTheta_in(j,:) + (0:nntheta:((nCluster-1)*nntheta));
        pick(opPick) = true;
        pick = reshape(pick,[nntheta,nCluster])';
        
        epirate(j,:) = CtC(j).frate(pick)';

        % orth
        pick = false(nCluster * nntheta,1);
        opPick = aiorthTheta_in(j,:) + (0:nntheta:((nCluster-1)*nntheta));
        pick(opPick) = true;
        pick = reshape(pick,[nntheta,nCluster])';
        
        eoirate(j,:) = CtC(j).frate(pick)';
    end
    clear opPick pick
    
    
    cc = contrastLevel;
    while cc <= 0
        cc = cc + 1;
    end
    
    shifted = 4;
    subplot(2,2,1); hold on
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eporate,2),std(eporate,1,2),'o:r');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoorate,2),std(eoorate,1,2),'o:k');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(epirate,2),std(epirate,1,2),'o:m');
    errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoirate,2),std(eoirate,1,2),'o:g');
    legend({'Prefout','Orthout','Prefin','Orthin'});
    
    pick = CnP(cc).PK > Cbr & CnP(contrastLevel).PK > thres;
    npick = sum(pick);
    if npick > 0
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'^-r');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'^-k');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'^-m');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'^-g');
    end
    ylim([0,inf]);xlim([0,110]);
    %legend('Prefered Orientation','Orthogonal Orientation')
    xlabel('Contrast Level %');
    ylabel('Firing Rate');
    title(['Exc ',num2str(npick/p.nv1e*100,'%3.1f'),'% evoked']);
    
    subplot(2,2,3); hold on
    
    pick = CnP(cc).SC(epick)>1.0;
    npickAll = sum(pick);
    if npickAll > 0
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'o:r');
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'o:k');
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'o:m');
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'o:g');
    end
    
    pick = CnP(cc).SC>1.0 & CnP(cc).PK > Cbr & CnP(contrastLevel).PK > thres;
    npick = sum(pick);
    if npick > 0
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'^-r');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'^-k');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'^-m');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'^-g');
    end
    
    ylim([0,inf]);xlim([0,110]);
    %legend('Prefered Orientation','Orthogonal Orientation')
    xlabel('Contrast Level %');
    ylabel('Firing Rate');
    title(['Simple Exc at ', conLabel{cc},' Contrast ', num2str(npick/npickAll*100,'%3.1f'),'% evoked']);
    
    subplot(2,2,4); hold on
    pick = CnP(cc).SC<1.0;
    npickAll = sum(pick);
    if npickAll > 0
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'o:r');
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'o:k');
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'o:m');
        errorbar(12.5.*2.^((1:contrastLevel)-1)+shifted,mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'o:g');
    end
    
    pick = CnP(cc).SC<1.0 & CnP(cc).PK > Cbr & CnP(contrastLevel).PK > thres;
    npick = sum(pick);
    disp([num2str(npick),' active complex cells at ',conLabel{cc},' contrast ', num2str(npickAll),' non-active complex cells']);
    if npick > 0
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eporate(:,pick),2),std(eporate(:,pick),1,2),'^-r');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoorate(:,pick),2),std(eoorate(:,pick),1,2),'^-k');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(epirate(:,pick),2),std(epirate(:,pick),1,2),'^-m');
        errorbar(12.5.*2.^((1:contrastLevel)-1),mean(eoirate(:,pick),2),std(eoirate(:,pick),1,2),'^-g');
    end
    
    ylim([0,inf]);xlim([0,110]);
    %legend('Prefered Orientation','Orthogonal Orientation')
    xlabel('Contrast Level %');
    ylabel('Firing Rate');
    title(['Complex Exc at ', conLabel{cc},' Contrast ', num2str(npick/npickAll*100,'%3.1f'),'% evoked']);
    
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hPrefOrth,[outputfdr,'/','GainCurve-Cluster-',theme,'.',format],printDriver,dpi);
    end
    
    clear eporate eoorate 
    if pOSI
        meanOSI = zeros(highestLGN-lowestLGN+1,contrastLevel)-1;
    %    stdOSI = meanOSI;
        aeOSI = zeros(contrastLevel,nCluster);
        aeOSInB = aeOSI;
        for i = 1:contrastLevel
            for j = 1:nCluster
                nominator = CtC(i).frate(j,aiprefTheta(i,j))-CtC(i).frate(j,aiorthTheta(i,j));
                denorm = CtC(i).frate(j,aiprefTheta(i,j))+CtC(i).frate(j,aiorthTheta(i,j));
                    aeOSI(i,j) = nominator/denorm;
                    aeOSInB(i,j) = nominator/(denorm -2*Cbr(j));
            end
        end
        aeOSI(isnan(aeOSI)) = 0;
        aeOSInB(isnan(aeOSInB)) = 0;
    
        hOSIpair = figure;
        subplot(1,2,1); hold on;
        p1 =contrastLevel-2; p2 = contrastLevel;
        while p1<=0
            p1 = p1+1;
        end
        %suptitle(theme);
        for i=lgnmin_e:lgnmax_e
            pick= CnLGN==i & CnP(p1).PK> Cbr & CnP(p2).PK> Cbr & CnP(contrastLevel).PK>thres;
            if ~isempty(pick)
                plot(aeOSI(p1,pick),aeOSI(p2,pick),'o','Color',[0.1+0.9*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
            end
            for j = 1:contrastLevel
                pick= CnLGN==i & CnP(j).PK> Cbr & CnP(contrastLevel).PK>thres;
                if ~isempty(pick)
                    meanOSI(i-lowestLGN+1,j) = mean(aeOSI(j,pick));
                end
            end
        end
        plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
        
        xlabel(['OSI at ',conLabel{p2},' Contrast']);
        ylabel(['OSI at ',conLabel{p1},' Contrast']);
    
        title('Excitatory');
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(hOSIpair,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_OSI-Cluster-',theme,'.',format],printDriver,dpi);
        end
        %% OSI over # LGN
        hOSIovernLGN_OSIvsCV = figure; 
    
        subplot(2,2,1); hold on;
        x = (-0.2:(highestLGN-lowestLGN-0.2));
        for i=1:contrastLevel
            %x = (-0.2:(highestLGN-lowestLGN-0.2))+0.1*(i-1);
            pick = meanOSI(:,i)>=0;
            plot(x(pick),meanOSI(pick,i),'-*','Color',[0.1+0.9*i/contrastLevel,0,0]);
        end
        xlabel('# LGN input'); ylabel('mean OSI');
        title('active Exc OSI');
        xlim([lgnmin_e, lgnmax_e]);
        ylim([0,inf]);
        yy = ylim();
        ylim([yy(1),yy(2)*1.2]);
        legend(conLabel);
    
        subplot(1,2,2); hold on;
        pp = 3;
        for i=lgnmin_e:lgnmax_e
            pick= CnLGN==i & CnP(pp).PK> Cbr & CnP(contrastLevel).PK>thres;
            plot(aeOSI(3,pick),1-CnP(pp).CV(pick),'o','Color',[0.1+0.9*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
        end
    
        plot(0:0.1:1,0:0.1:1,'-.k','LineWidth',2);
        xlabel('OSI');
        ylabel('1-CV');
        title(['1-CV vs. OSI at',num2str(12.5*2^(pp-1),'%3.1f'),'% Contrast']);
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            print(hOSIovernLGN_OSIvsCV,[outputfdr,'/','OSIovernLGN_OSIvsCV-Cluster-',theme,'.',format],printDriver,dpi);
        end
        clear aeOSI aeOSInB
    end
    %% simple complex
    hSimCom = figure;
    meanSC = zeros(highestLGN-lowestLGN+1,contrastLevel)-1;
    
    p1 =contrastLevel-2; p2 = contrastLevel;
    while p1<=0
        p1 = p1+1;
    end
    hold on
    for i=lgnmin_e:lgnmax_e
        pick= CnP(contrastLevel).ei>0.5 &CnLGN==i & CnP(p1).PK> Cbr & CnP(p2).PK> Cbr & CnP(contrastLevel).PK>thres;
        if ~isempty(pick)
            plot(CnP(p1).SC(pick),CnP(p2).SC(pick),'o','Color',[0.1+0.9*(i-lgnmin_e)/(lgnmax_e-lgnmin_e),0,0]);
        end
        for j = 1:contrastLevel
            pick= CnLGN==i & CnP(j).PK> Cbr & CnP(contrastLevel).PK>thres;
            if ~isempty(pick)
                meanSC(i-lowestLGN+1,j) = mean(CnP(j).SC(pick));
            end
        end
    end
    x = 0:0.1:2;
    nx = length(x);
    plot(x,linspace(0,2,nx),'-.k','LineWidth',2);
    plot(x,ones(nx,1),'-.k','LineWidth',2);
    plot(ones(nx,1),x,'-.k','LineWidth',2);
    
    xlabel(['F1/F0 at ',conLabel{p2},' Contrast']);
    ylabel(['F1/F0 at ',conLabel{p1},' Contrast']);
    
    title('Excitatory');
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hSimCom,[outputfdr,'/','C',num2str(p1),'vsC',num2str(p2),'_F1overF2-Cluster-',theme,'.',format],printDriver,dpi);
    end
    clear meanSC
    %% orientations
    hOriExc = figure;
    
    thetaRange = 180/ntheta*((0:ntheta)-0.5);
    for j = 1:contrastLevel
        pick = CnP(j).PK> Cbr & CnP(contrastLevel).PK>thres;
        oriRate = zeros(ntheta,1);
        po1d = (CnP(j).iPrA(pick)-1)*dtheta;
        subplot(contrastLevel,2,(j-1)*2+1);
        [binCount, ind] = histc(po1d,thetaRange);
        binCount = binCount(1:ntheta);
        fired_PK = CnP(j).PK(pick);
        for i=1:ntheta
            thetaPick = (ind == i);
            if sum(thetaPick)>0
                oriRate(i) = mean(fired_PK(thetaPick));
            end
        end
        plot((0:(ntheta-1))*180/ntheta,oriRate,'LineWidth',2);
        xlim([0,180-180/ntheta]);
        xlabel('pref angle (^{o})');
        ylabel('Hz');
    
        subplot(contrastLevel,2,(j-1)*2+2);
        plot((0:(ntheta-1))*180/ntheta,binCount,'LineWidth',2);
        xlim([0,180-180/ntheta]);
        xlabel('pref angle (^{o})');
        ylabel('# Exc Cell');
        title([num2str(sum(pick)/p.nv1e*100,'%2.1f'),'% Exc evoked']);
    end
    
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hOriExc,[outputfdr,'/','ExcOrientation-Cluster-',theme,'.',format],printDriver,dpi);
    end
    hOriCmp = figure; hold on;
    p1 = 3; p2 = 4;
    shifts = zeros(lgnmax_e-lgnmin_e+1,1);
    subplot(1,2,1); hold on

    pick= CnP(contrastLevel).PK> Cbr & CnP(contrastLevel).PK>thres;
    [y,x]=hist((aiprefTheta(contrastLevel,:)-1)*dtheta,-dtheta:dtheta:180-dtheta);
    bar(x,y);
    [y,x]=hist((aiprefTheta(contrastLevel,pick)-1)*dtheta,-dtheta:dtheta:180-dtheta);
    hBar = bar(x,y);
    set(hBar,'FaceColor','r');
    legend({'Inh','Exc'});
    xlabel('Preferred Orientation at the Highest Contrast (degree)');    
    xlim([-dtheta,180+dtheta]);
    for i=lgnmin_e:lgnmax_e
        pick= CnLGN==i & CnP(contrastLevel).PK> Cbr & CnP(contrastLevel).PK>thres;
        if ~isempty(pick)
            shifts(i-lgnmin_e+1) = sum(abs(aiprefTheta(p2,pick)-aiprefTheta(p1,pick)))/sum(pick)*dtheta;
        end
    end
    
    
    subplot(1,2,2);
    bar((lgnmin_e+0.5):(lgnmax_e+0.5),shifts);
    xlabel('# LGN input');
    ylabel('mean shift in preferred orientation (^{o})');
    title('Exc only');
    xlim([lgnmin_e,inf]);
    if ~isempty(format)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        print(hOriCmp,[outputfdr,'/','OriCmp-Cluster-',theme,'.',format],printDriver,dpi);
    end
    end
    disp('population statistics done');
    toc;
end
%%
function [X, Y, Z] = multiLGNspatialKernel(posx,posy,x,y,n,s,scale)
    Aa = 14.88 * (180/pi)^2;
    Ab = Aa*0.97;
    rc = 5.61; %deg
    rs = 16.98*rc/5.61; %deg
    siga = rc/sqrt(2);
    sigb = rs/sqrt(2);
    Z = zeros(length(y),length(x));
    [X, Y] = meshgrid(x,y);
    for i = 1:n
        Z = Z + s(i)*(Aa/(siga^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/siga^2/2)...
            - Ab/(sigb^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/sigb^2/2));
    end
    Z = Z/scale;
end
