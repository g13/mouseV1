function [p, etheta, itheta, sp, nLGN, nSubLGN, v1Map, pLGN, lgnStrength,v1pos] = lgn2v1Map_beta(p)
    global dLGNpos LGNpos pick rfx rfy ttt rlgn rho gauss se si sbounde sboundi cbounde cboundi dlgn uniform profile printDriver dpi pPosition uniformStrength FWHM HWHM Yidx
    FontSize = 14;
    set(0,'DefaultAxesFontSize',FontSize);
    pPosition = [18, 180, 1200, 900];
    if ~isempty(p.h)
        if strcmp(p.h,'psc2')
            printDriver = ['-de',p.h];
            p.h = 'eps';
        else
            printDriver = ['-d',p.h];
        end
        dpi = '-r150';
    end
	% position in rads
    rlgn = p.rlgn;
    pick = false(p.nlgn,1);
    LGNpos = zeros(p.nlgn,2);
    % coordinate x / column
	LGNpos(:,1) = reshape( ones(p.lgny,1) * linspace(p.lgn_azimuth(1),p.lgn_azimuth(2),p.lgnx), [p.nlgn,1]);
	% coordinate y / row
	LGNpos(:,2) = reshape( (linspace(p.lgn_altitude(1),p.lgn_altitude(2),p.lgny))' * ones(1,p.lgnx), [p.nlgn,1]);
    dlgnx = (p.lgn_azimuth(2) - p.lgn_azimuth(1))/(p.lgnx-1);
    dlgny = (p.lgn_altitude(2) - p.lgn_altitude(1))/(p.lgny-1);
    disp(['x-axis lgn distance: ', num2str(dlgnx*180/pi,'%3.1f')]);
    disp(['y-axis lgn distance: ', num2str(dlgny*180/pi,'%3.1f')]);
    disp(['lgn radius: ', num2str(rlgn*180/pi,'%3.1f')]);
    
    cbounde = p.cbounde;
    cboundi = p.cboundi;
    sbounde = p.sbounde;
    sboundi = p.sboundi;
    se = p.se;
    si = p.si;
    FWHM = 2*sqrt(2*log(2));
    HWHM = FWHM/2;
    
%     rho = @(slope,a,b) a.^2.*b.^2.*( (1+slope.^2)./(b.^2+a.^2.*slope.^2) );
    rho = @(theta,a,b) (a.*b).^2./((b.*cos(theta)).^2+(a*sin(theta)).^2);
    % draw lgnrf (ellipse and circle)
    rfx = @(x0,a,b,theta,t,aratio,bratio) x0 + (aratio.*a.*cos(t).*cos(theta)-bratio.*b.*sin(t).*sin(theta));
    rfy = @(y0,a,b,theta,t,aratio,bratio) y0 + (aratio.*a.*cos(t).*sin(theta)+bratio.*b.*sin(t).*cos(theta));
    ttt = 0:0.1:2*pi;
    Yidx = p.YoudensIndex;
%     seed = p.seed;
    dlgn = p.dlgn;
    profile = p.profile;
    if p.seed >=0
        if p.seed ==0
            rng('shuffle');
        else
            rng(p.seed);
        end
        if profile == 'u' || profile == 'g'
%             gauss = @(x,y,a,b) dlgnx*dlgny*bound^2*dlgn/2*...
%                     exp(-0.5*(x.^2./a.^2+y.^2./b.^2));
            gauss = @(x,y,a,b) dlgnx*dlgny*dlgn/(2*pi*a*b)*...
                    exp(-0.5*(x.^2./a.^2+y.^2./b.^2));
            uniform = dlgn*dlgnx*dlgny;
            if uniform > 1
                uniform = 1;
                disp('prescribed density require denser lgn distribution');
            end
%             dlgn = uniform/(dlgnx*dlgny);
            if profile == 'u'
%                 p.average = 1; % no repeat
                disp(['uniform prob: ',num2str(uniform)]);
            end
        else
            disp(['no profile is defined for ''',p.profile,'''']);
            return
        end
    else
        disp('seed must be positive');
        return
    end
    uniformStrength = p.uniformStrength;
    if p.lgn_offset ~= 0
        xoffset = randn([p.nlgn,1]) * p.lgn_offset * (LGNpos(p.lgny+1,1)-LGNpos(1,1));
        yoffset = randn([p.nlgn,1]) * p.lgn_offset * (LGNpos(2,2)-LGNpos(1,2));
    end
    
    LGNpos(:,1) = LGNpos(:,1) + xoffset;
    LGNpos(:,1) = LGNpos(:,1) + yoffset;
    dLGNpos = LGNpos.*180./pi;
    h = figure;
	subplot(2,2,1);
    hold on
    plot( dLGNpos(:,1), dLGNpos(:,2),'x','MarkerSize',2.0);
    xlabel('degree'); ylabel('degree');

    dev1y = (p.v1_altitude(2)-p.v1_altitude(1))/p.ev1y; % for correction upon periodical boundary
    axis_ev1y = linspace(p.v1_altitude(1)+0.5*dev1y,p.v1_altitude(2)-0.5*dev1y,p.ev1y);
 
    dev1x = (p.v1_azimuth(2)-p.v1_azimuth(1))/p.ev1x;
    axis_ev1x = linspace(p.v1_azimuth(1)+0.5*dev1x,p.v1_azimuth(2)-0.5*dev1x,p.ev1x);
    
    ev1posx = ones(p.ev1y,1) * axis_ev1x;
    ev1posy = axis_ev1y' * ones(1,p.ev1x);
    
    div1y = (p.v1_altitude(2)-p.v1_altitude(1))/p.iv1y;
    axis_iv1y = linspace(p.v1_altitude(1)+0.5*div1y,p.v1_altitude(2)-0.5*div1y,p.iv1y);
    
    div1x = (p.v1_azimuth(2)-p.v1_azimuth(1))/p.iv1x;
    axis_iv1x = linspace(p.v1_azimuth(1)+0.5*div1x,p.v1_azimuth(2)-0.5*div1x,p.iv1x);
    
    iv1posx = ones(p.iv1y,1) * axis_iv1x;
    iv1posy = axis_iv1y' * ones(1,p.iv1x);
    
    title({['\Delta^{E}_{x} = ',num2str(dev1x*180/pi,'%3.1f'),', \Delta^{E}_{y} = ',num2str(dev1y*180/pi,'%3.1f')];...
        ['\Delta^{I}_{x} = ',num2str(div1x*180/pi,'%3.1f'),', \Delta^{I}_{y} = ',num2str(div1y*180/pi,'%3.1f')];...
        ['\Delta^{LGN}_{x} = ',num2str(dlgnx*180/pi,'%3.1f'),', \Delta^{LGN}_{y} = ',num2str(dlgny*180/pi,'%3.1f')]});
    axis equal
    % assign orientation
    if ~isempty(p.pinwheel)
        etheta = zeros(p.ev1y,p.ev1x);
        itheta = zeros(p.iv1y,p.iv1x);
        dxcenter = diff(p.v1_azimuth)/(p.pinwheel(2));
        dycenter = diff(p.v1_altitude)/(p.pinwheel(1));
        cev1x = p.ev1x/p.pinwheel(2);
        civ1x = p.iv1x/p.pinwheel(2);
        cev1y = p.ev1y/p.pinwheel(1);
        civ1y = p.iv1y/p.pinwheel(1);
        for ix=1:p.pinwheel(1)
            for iy = 1:p.pinwheel(2)
                pincenter(1) = p.v1_azimuth(1) + (ix-1/2)*dxcenter;
                pincenter(2) = p.v1_altitude(1) + (iy-1/2)*dycenter;
                for i = 1:cev1x
                    for j = 1:cev1y
                        ipos = min(cev1x*(ix-1) + i,p.ev1x);
                        jpos = min(cev1y*(iy-1) + j,p.ev1y);
                        etheta(jpos,ipos) = atan((ev1posy(jpos,ipos)-pincenter(2))/(ev1posx(jpos,ipos)-pincenter(1)));
                        if isnan(etheta(jpos,ipos))
                            etheta(jpos,ipos) = pi/2;
                        end
                    end
                end
                for i = 1:civ1x
                    for j = 1:civ1y
                        ipos = min(civ1x*(ix-1) + i,p.iv1x);
                        jpos = min(civ1y*(iy-1) + j,p.iv1y);
                        itheta(jpos,ipos) = atan((iv1posy(jpos,ipos)-pincenter(2))/(iv1posx(jpos,ipos)-pincenter(1)));
                        if isnan(itheta(jpos,ipos))
                            itheta(jpos,ipos) = pi/2;
                        end
                    end
                end
            end
        end
        etheta(etheta<0) = etheta(etheta<0) + pi;
        itheta(itheta<0) = itheta(itheta<0) + pi;
    else
        etheta = pi*rand([p.ev1y,p.ev1x]);
        itheta = pi*rand([p.iv1y,p.iv1x]);
    end
    ev1pos = [reshape(ev1posx,[p.nv1e,1]),reshape(ev1posy,[p.nv1e,1])];
    iv1pos = [reshape(iv1posx,[p.nv1i,1]),reshape(iv1posy,[p.nv1i,1])];
    v1pos = [ev1pos; iv1pos];
    [v1Map, nLGN, nSubLGN, sp, v1pos, lgnStrength, ORF, itermax, etheta, itheta, iCRF] = gRF(p,etheta,itheta,v1pos);
    disp(['maximum iterations: ',num2str(max(itermax))]);
    p.sp = sp;
    sp = sp/pi*180;
    p.ORF = ORF;
    p.CRF = iCRF;
    % specifies ORF, SRF and CRF Exc
    pick2e = [p.eSubregion == 2; false(p.nv1i,1)];
    p.typeE(ORF & pick2e) = 1;
    p.typeE(~ORF & pick2e) = 2; 
    p.typeE(iCRF & pick2e) = 3;
    assert(sum(ORF-iCRF<0)==0);
    clear pick2e;
    %% plot
    figure(h);
    subplot(2,2,1);
    hold on
    % draw colormap of orientation RF
    if p.drawOrientation
        thet = reshape(etheta,[p.nv1e,1]);
        thet(thet>=pi) = thet(thet>=pi)-pi;
        thet(thet<0) = thet(thet<0)+pi;
        hsv = zeros(p.nv1e,3);
        hsv(:,1) = thet./pi;
        hsv(:,2) = 0.7;
        hsv(:,3) = 0.9;
        color = hsv2rgb(hsv);
        for i=1:p.nv1e
            plot(v1pos(i,1)/pi*180, v1pos(i,2)/pi*180 ,'.','Color',color(i,:),'MarkerSize',8.0); 
        end

        thet = reshape(itheta,[p.nv1i,1]);
        thet(thet>=pi) = thet(thet>=pi)-pi;
        thet(thet<0) = thet(thet<0)+pi;
        hsv = zeros(p.nv1i,3);
        hsv(:,1) = thet./pi;
        hsv(:,2) = 0.7;
        hsv(:,3) = 0.9;
        color = hsv2rgb(hsv);
        for i=1:p.nv1i
            plot(v1pos(p.nv1e+i,1)/pi*180, v1pos(p.nv1e+i,2)/pi*180,'.','Color',color(i,:),'MarkerSize',8.0); 
        end
    else
         plot(v1pos(1:p.nv1e,1)/pi*180,v1pos(1:p.nv1e,2)/pi*180,'.r','MarkerSize',8.0);
         plot(v1pos(p.nv1e+1:end,1)/pi*180,v1pos(p.nv1e+1:end,2)/pi*180,'.b','MarkerSize',8.0);
    end

    % nLGN distribution
    id = (1:p.nv1)';
    subplot(2,2,2);
    hold on;
    centers = min(nLGN(:,1)):5:max(nLGN(:,1));
    hist(nLGN(:,1),centers);
    hhist = findobj(gca,'Type','patch');
    hhist.FaceColor = [0 0.5 0.0];
    hhist.EdgeColor = 'w';
    subregion = [p.eSubregion; p.iSubregion];
    mS = max(subregion);
    binCounts = zeros(mS+1,length(centers));
    for k = 1:mS
         binCounts(k,:) = histc(nLGN(subregion==k & id<=p.nv1e,1),centers);
    end
    binCounts(mS+1,:) = histc(nLGN(id>p.nv1e,1),centers);
    plot(centers,binCounts,'LineWidth',2.0);
    xlim([0,inf]);
    xlabel('# LGN input');
    legend(['total';cellstr(num2str((1:mS)'));'inh']);
    meanLGNe = mean(nLGN(1:p.nv1e,1));
    meanLGNi = mean(nLGN(p.nv1e+(1:p.nv1i),1));
    title({['mean(nLGN_E) = ',num2str(meanLGNe,'%3.1f')],['mean(nLGN_I) = ',num2str(meanLGNi,'%3.1f')]});
    
    % spatial frequency
    subplot(2,2,3);
    hold on
    
    centers = linspace(min(sp),max(sp),10);
    hist(sp,centers);
    hhist = findobj(gca,'Type','patch');
    hhist.FaceColor = [0 0.5 0.0];
    hhist.EdgeColor = 'w';
    
    binCounts = zeros(mS+1,length(centers));
    for k = 1:mS
        binCounts(k,:) = histc(sp(subregion==k & id<=p.nv1e),centers);
    end
    binCounts(mS+1,:) = histc(sp(id>p.nv1e),centers);
    plot(centers,binCounts,'LineWidth',2.0);
    xlim([0,inf]);
    
    ORFcount = hist(sp([p.typeE == 1; false(p.nv1i,1)]),centers);
    SRFcount = hist(sp([p.typeE == 2; false(p.nv1i,1)]),centers);
    plot(centers,[ORFcount; SRFcount],':','LineWidth',2.0);
    xlim([0,inf]);
    
    legend(['total';cellstr(num2str((1:mS)'));'inh';'ORF';'SRF']);
    xlabel('Spatial Period (degree, 1 by FWHM, >1 by peak distance)'); ylabel('count');
    n2sub = sum(p.eSubregion == 2);
    disp([num2str(sum(p.typeE==1)/n2sub*100,'%3.1f'),'% ORF vs ',...
        num2str(sum(p.typeE==2)/n2sub*100,'%3.1f'), '% SRF']);
    title([num2str(sum(p.typeE==1)/n2sub*100,'%3.1f'),'% ORF vs ',...
        num2str(sum(p.typeE==2)/n2sub*100,'%3.1f'), '% SRF']);
    
%     [~, ind] = sort(binCounts,'descend');
%     f4 = ind(1:4);
%     opt_sp = binCounts(f4)*centers(f4)'/sum(binCounts(f4));
%     opt_gk = 1/opt_sp;
% 	disp(['optimal spatial period: ',num2str(opt_sp)]);

    subplot(2,2,4);
    hold on
    centers = min(nLGN(:,4)):5:max(nLGN(:,4));
    hist(nLGN(:,4),centers);
    hhist = findobj(gca,'Type','patch');
    hhist.FaceColor = [0 0.5 0.0];
    hhist.EdgeColor = 'w';
    subregion = [p.eSubregion; p.iSubregion];
    mS = max(subregion);
    binCounts = zeros(mS+1,length(centers));
    for k = 1:mS
        binCounts(k,:) = histc(nLGN(subregion==k & id<=p.nv1e,4),centers);   
    end
    binCounts(mS+1,:) = histc(nLGN(id>p.nv1e,4),centers);
    plot(centers,binCounts,'LineWidth',2.0);
    xlim([0,inf]);
    xlabel('# LGN input');
    legend(['total';cellstr(num2str((1:mS)'));'inh']);
    meanLGNe = mean(nLGN(1:p.nv1e,4));
    meanLGNi = mean(nLGN(p.nv1e+(1:p.nv1i),4));
    title({['Non-overlap mean(nLGN_E) = ',num2str(meanLGNe,'%3.1f')],['mean(nLGN_I) = ',num2str(meanLGNi,'%3.1f')]});
%%%%%%%% used to print estimated frequency based on peak distance, rather
%%%%%%%% should be estimated by subregion size
%     subplot(2,2,4);
%     hold on
%     centers = 0:12;
%     factor = (2^(1/3));
%     tickLabel = [1,2,4,8,16]';
%     tickPosition = centers(1:3:13);
%     pos = 0.01*factor.^(centers);
% %     cutoff = norm(dlgnx,dlgny);
%     gk = 1./sp;
%     gk = gk(~isnan(gk));
%     binCount = hist(gk(~ORF),pos);
%     [f3,i3] = sort(binCount,'descend');
%     f3 = f3(1:3); i3 = i3(1:3);
%     p.opt_gk = sum(((pos(i3)) .* f3))/sum(f3);
%     disp(['estimated optimal sp freq: ',num2str(p.opt_gk*180/pi)]);
%     bar(centers,binCount);
%     binCounts = zeros(mS+1,length(centers));
%     for k=1:mS
%         binCounts(k,:) = hist(gk(~ORF & subregion == k & id<=p.nv1e),pos);
%     end
%     binCounts(mS+1,:) = hist(gk(~ORF & id > p.nv1e),pos);
%     plot(centers,binCounts,'LineWidth',2.0);
%     xlim([0,centers(end)+1]);
%     legend(['total';cellstr(num2str((1:mS)'));'inh']);
%     set(gca,'XTick',tickPosition);
%     set(gca,'XTickLabel',cellstr(num2str(tickLabel)));
% %     hist(gk,0:15);
% %     xlabel('Spatial Frequency (0Hz -> ORF)');
%     title({['optimal spatial freq',num2str(p.opt_gk),'cpd'],...
%            '(estimated over the top 3 bars, ORF excluded)'});
%     xlabel('Spatial Frequency x 10^{-2} (cpd, by peak distance)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(p.h)
        set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
        if strcmp(p.h,'fig')
            savefig(h,[p.fdr,'/lgn2v1.',p.h]);
        else
            print(h,[p.fdr,'/lgn2v1.',p.h],printDriver,dpi);
        end
    end
    %% output to lgnmat.out for simulation
    fid = fopen(['lgn2v1map/lgnmap','.out-',p.name],'w+');
    lofexplanation = 9;
    fprintf(fid,'*********************** \n');
    if p.pick
        lofLGN = sum(pick);
    else
        lofLGN = p.nlgn;
    end
    fprintf(fid,'%i LGNs'' x,y position from line %i to %i.\n',lofLGN, lofexplanation + 1, lofexplanation + lofLGN);
    fprintf(fid,'----------------------- \n');
    fprintf(fid,'V1 cells'' info start from line %i with the following format:\n',lofexplanation + lofLGN + 1);
    fprintf(fid,'#sub regions, prescribed preferred drifting orientation (orthogonal to the grating orientation)\n');
    fprintf(fid,'#LGN in the ith sub region\n');
    fprintf(fid,'id of the jth LGN in ith sub region, strength of the jth LGN\n');
    fprintf(fid,'loop over all sub regions and then follows the next V1 cell...\n');
    fprintf(fid,'*********************** \n');
    if p.pick
		p.npick=sum(pick);
        disp(['Number of lgn pick: ',num2str(p.npick),', ',num2str(100*p.npick/p.nlgn),'%']);
        fprintf(fid,'%G %G \n',LGNpos(pick,:));
		p.ipick = zeros(p.nlgn,1);
		p.ipick(pick) = 1:p.npick;
        assert(iscolumn(p.ipick));
    else
        fprintf(fid,'%G %G \n',LGNpos');
    end
    thet = [reshape(etheta,[p.nv1e,1]);reshape(itheta,[p.nv1i,1])];
    nSubregion = [p.eSubregion; p.iSubregion];
    
    p.EmeanLGN = mean(nLGN(1:p.nv1e,1));
    p.ImeanLGN = mean(nLGN(p.nv1e+(1:p.nv1i),1));
%%% homogeneous LGN input
%    lgninput = zeros(p.nv1,1);
%%%
    if p.pick
        for i=1:p.nv1
            fprintf(fid,'%i %G \n', nSubregion(i), thet(i));
            for j = 1:nSubregion(i)
                fprintf(fid,'%i \n',nSubLGN{i}(j));
%%% homogeneous LGN input
%    if i<=p.nv1e
%        lgnStrength{i,j} = lgnStrength{i,j}.*p.EmeanLGN/nLGN(i,1);
%        lgninput(i) = lgninput(i) + sum(lgnStrength{i,j});
%    else
%        lgnStrength{i,j} = lgnStrength{i,j}.*p.ImeanLGN/nLGN(i,1);
%        lgninput(i) = lgninput(i) + sum(lgnStrength{i,j});
%    end
%%%%
                if nSubLGN{i}(j) < 0
                    fprintf(fid,'%i %G \n',[-p.ipick(-v1Map{i,j})';lgnStrength{i,j}']);
                    assert(iscolumn(-p.ipick(-v1Map{i,j})));
                    assert(iscolumn(lgnStrength{i,j}));
                else
                    fprintf(fid,'%i %G \n',[p.ipick(v1Map{i,j})';lgnStrength{i,j}']);
                    assert(iscolumn(p.ipick(v1Map{i,j})));
                    assert(iscolumn(lgnStrength{i,j}));
                end
                
            end
        end
    else
        for i=1:p.nv1
            fprintf(fid,'%i %G \n', nSubregion(i), thet(i));
            for j = 1:nSubregion(i)
%%% homogeneous LGN input
%    if i<=p.nv1e
%        lgnStrength{i,j} = lgnStrength{i,j}.*p.EmeanLGN/nLGN(i,1);
%        lgninput(i) = lgninput(i) + sum(lgnStrength{i,j});
%    else
%        lgnStrength{i,j} = lgnStrength{i,j}.*p.ImeanLGN/nLGN(i,1);
%        lgninput(i) = lgninput(i) + sum(lgnStrength{i,j});
%    end
%%%
                fprintf(fid,'%i \n',nSubLGN{i}(j));
                fprintf(fid,'%i %G \n',[v1Map{i,j}; lgnStrength{i,j}']);
                assert(isrow(v1Map{i,j}));
                assert(iscolumn(lgnStrength{i,j}));
            end
        end
    end
    pLGN = LGNpos;
    p.lgnmin = min(nLGN(:,1));
    fclose(fid);
    disp(['# Complete ORF = ',num2str(sum(iCRF))]);
%%% homogeneous LGN input
%    figure
%    hist(lgninput,10);
%    title([num2str(std(lgninput(1:p.nv1e))),',',num2str(std(lgninput(p.nv1e+(1:p.nv1i))))]);
%%%
end

function [v1Map, nLGN, nSubLGN, sp, v1pos, lgnStrength, iORF, itermax, etheta, itheta, iCRF] = gRF(p,etheta,itheta,v1pos)
    global pick cbounde cboundi se si sbounde sboundi
	
    etheta = reshape(etheta,[p.nv1e,1]);
    maxSubgregion = max([p.eSubregion;p.iSubregion]);
    v1Map = cell(p.nv1,maxSubgregion);
    lgnStrength = cell(p.nv1,maxSubgregion);
    sp = zeros(p.nv1,1);
    iORF = false(p.nv1,1);
    iCRF = false(p.nv1,1);
    nLGN = zeros(p.nv1,4);
    nSubLGN = cell(p.nv1,1);
    itermax = zeros(p.nv1,1);
    example = 0;
    if p.example > 16
        set(0,'DefaultFigureVisible','off');
    end
    for i = 1:p.nv1e
%         if mod(i,100) == 1
%             disp(i);
%         end
        if p.eSubregion(i) == 0
            v1Map(i) = {};
            sp(i) = inf;
        else
            eex = ceil(p.example*p.nv1e/p.nv1);
            if example < eex && rand < eex/p.nv1e
                example = example+1;
                ex = example;
            else
                ex = 0;
            end
            theta = etheta(i) + randn(p.eSubregion(i),1)*p.sigmaTheta;
            
            q = roll([num2str(i),' ',p.Etypes{p.typeE(i)}],...
                    p.eSubregion(i),p.esigma(i,1:p.eSubregion(i)),...
                    p.eAspectRatio(i,1:p.eSubregion(i)), ...
                    p.ePeakDistance(i,1:max(1,p.eSubregion(i)-1)),etheta(i),...
                    theta,p.nlgn,v1pos(i,:),ex,p.pOn,p.h,p.name,p.average,cbounde,se,sbounde,true);
            if length(q.v1Map{1}) == length(q.v1Map{2})
                if sum(q.v1Map{1}+q.v1Map{2}) == 0
                    iCRF(i) = true;
                end
            end
            v1Map(i,1:p.eSubregion(i)) = q.v1Map; 
            lgnStrength(i,1:p.eSubregion(i)) = q.strength;
            v1pos(i,:) = q.center;
            sp(i) = q.sp;
            pick(q.pick) = true;
            nLGN(i,:) = [q.nLGN,q.nLGNon,q.nLGNoff,sum(q.pick)];
            nSubLGN{i} = q.isub;
            iORF(i) = q.iORF;
            itermax(i) = q.itermax;
            etheta(i) = q.theta;
        end
    end
    etheta(etheta<0) = etheta(etheta<0) + pi;
    etheta(etheta>=pi) = etheta(etheta>=pi) - pi;
    
    itheta = reshape(itheta,[p.nv1i,1]);
    example = 0;
    for i = 1:p.nv1i
        if p.iSubregion(i) == 0
            v1Map(p.nv1e+i) = {};
            sp(i) = inf;
        else 
            iex = ceil(p.example*p.nv1i/p.nv1);
            if example < iex && rand < iex/p.nv1i
                example = example+1;
                ex = example;
            else
                ex = 0;
            end
            theta = itheta(i) + randn(p.iSubregion(i),1)*p.sigmaTheta;
            q = roll([num2str(i+p.nv1e),' ',p.Itypes{p.typeI(i)}],...
                    p.iSubregion(i),p.isigma(i,1:p.iSubregion(i)),...
                    p.iAspectRatio(i,1:p.iSubregion(i)), ...
                    p.iPeakDistance(i,1:max(1,p.iSubregion(i)-1)),itheta(i),...
                    theta,p.nlgn,v1pos(p.nv1e+i,:),ex,p.pOn,p.h,p.name,p.average,cboundi,si,sboundi,false);

            v1Map(p.nv1e+i,1:p.iSubregion(i)) = q.v1Map;
            lgnStrength(p.nv1e+i,1:p.iSubregion(i)) = q.strength;
            v1pos(p.nv1e+i,:) = q.center;
            sp(p.nv1e+i) = q.sp;
            pick(q.pick) = true;
            nLGN(p.nv1e+i,:) = [q.nLGN,q.nLGNon,q.nLGNoff,sum(q.pick)];
			nSubLGN{p.nv1e+i} = q.isub;
            iORF(p.nv1e+i) = q.iORF;
            itermax(p.nv1e+i) = q.itermax;
            itheta(i) = q.theta;
        end
    end
    itheta(itheta<0) = itheta(itheta<0) + pi;
    itheta(itheta>=pi) = itheta(itheta>=pi) - pi;
    if p.example > 16
        set(0,'DefaultFigureVisible','on');
    end
end
function q = roll(type,subregion,sigma,aspectRatio,peakDistance,prescribedTheta,theta,nlgn,center,example,pOn,format,suffix,nrep,cbound,s,sbound,ei)
    global LGNpos dLGNpos rfx rfy ttt rho rlgn gauss dlgn uniform profile printDriver dpi pPosition uniformStrength FWHM HWHM Yidx
    % theta is aligned with drifting direction, orthogonal to the orientation of the grating 
    gtheta = theta;
    pick = theta>=pi/2;
    % gtheta is the orientation of the grating
    gtheta(pick)= theta(pick) -pi/2;
    gtheta(~pick) = theta(~pick) + pi/2;
    cosgtheta = cos(gtheta);
    singtheta = sin(gtheta);
    opeak = zeros(subregion,2);
    peak = opeak;
%     slope = tan(theta);
    q.pick = false(nlgn,1);
    q.nLGNon = 0;
    q.nLGNoff = 0;
    q.nLGN = 0;
    q.iORF = false;
    q.sp = inf;
    
    q.v1Map = cell(subregion,1);
    q.strength = cell(subregion,1);
    q.isub = zeros(subregion,1);
    
    
    if subregion == 1
        opeak = center;
    else
        d = cumsum(peakDistance);
        opeak(1,1) = center(1) - d(end)/2*cos(prescribedTheta);
        opeak(1,2) = center(2) - d(end)/2*sin(prescribedTheta);
        opeak(2:subregion,1) = opeak(1,1) + d*cos(prescribedTheta);
        opeak(2:subregion,2) = opeak(1,2) + d*sin(prescribedTheta);
    end
    % note that sigmb is the semi major axis.
    bounda = cbound*ones(1,subregion);
    if ei
        boundb = bounda.*aspectRatio;
    else
        boundb = bounda.*aspectRatio;
%         boundb = bounda;
    end
%     boundThres = 8/(2*sqrt(2*log(2)))/180*pi;
%     bounda(sigma<boundThres) = bound/sqrt(2)/1.3;
    
    sigmb = sigma.*aspectRatio;
    prescribedNeff = dlgn*pi*sigma.*bounda.*sigmb.*boundb;
    nactual = ceil(prescribedNeff);

%     roll = rand(nlgn,subregion);
    q.itermax = 0;
    id = 1:nlgn;
    
    for k=1:subregion
        d2c2 = (LGNpos(:,2)-opeak(k,2)).^2 + (LGNpos(:,1)-opeak(k,1)).^2;
        ptheta = atan((LGNpos(:,2) - opeak(k,2))./(LGNpos(:,1)-opeak(k,1)));
        
        pick = ptheta < 0;
        ptheta(pick) = pi + ptheta(pick) - gtheta(k);
        ptheta(~pick) = ptheta(~pick) - gtheta(k);
        r2 = rho(ptheta,boundb(k)*sigmb(k),bounda(k)*sigma(k));
%         assert(~isnan(sum(r2)));
        valid = (d2c2 < r2);
        nvalid = sum(valid);
        idvalid = id(valid);
        switch profile
            case 'g'
                % canonical coordinates regard to the center 
                xcan = (LGNpos(valid,1)-opeak(k,1)).*cosgtheta(k) + (LGNpos(valid,2)-opeak(k,2)).*singtheta(k);
                ycan = -(LGNpos(valid,1)-opeak(k,1)).*singtheta(k) + (LGNpos(valid,2)-opeak(k,2)).*cosgtheta(k);
                
%                 xcan = (LGNpos(valid,1)-peak(k,1)).*sintheta(k) + (LGNpos(valid,2)-peak(k,2)).*costheta(k);
%                 ycan = -(LGNpos(valid,1)-peak(k,1)).*costheta(k) + (LGNpos(valid,2)-peak(k,2)).*sintheta(k);
                 
                cap = gauss(xcan,ycan,sigmb(k),sigma(k));
                cap = cap*ones(1,nrep);
            case 'u'
                cap = uniform;
        end
        % repeat dice roll
        if	nvalid > 0
            iter = 1;
            if nvalid <= nactual(k)
                nactual(k) = nvalid;
                picked = idvalid;
            else
                while true
                    roll = rand(nvalid,nrep);
                    result = mean(roll < cap ,2);
                    [~, rank] = sort(result,'descend');
    %                 if sum(result>0) >= neff(k)*0.9
                    if sum(result>0) >= nactual(k)
                        picked = idvalid(rank(1:nactual(k)));
                        break
                    end
                    iter = iter + 1;
                    nrep = nrep+nactual;
                end
            end
        else
            error('need denser LGN');
        end
        if iter > q.itermax;
            q.itermax = iter;
        end
        if nvalid > 0
			q.pick(picked) = true;
            q.isub(k) = nactual(k);
			q.v1Map{k} = picked;
        end
        q.strength{k} = zeros(nactual(k),1);
%         peak(k,:) = mean(LGNpos(q.v1Map{k},:),1);
        peak = opeak;

        if uniformStrength
            q.strength{k} = ones(q.isub(k),1);
        else 
            xcan = (LGNpos(picked,1)-peak(k,1)).*cosgtheta(k) + (LGNpos(picked,2)-peak(k,2)).*singtheta(k);
            ycan = -(LGNpos(picked,1)-peak(k,1)).*singtheta(k) + (LGNpos(picked,2)-peak(k,2)).*cosgtheta(k);
            q.strength{k} = gauss(xcan,ycan,sbound*sigmb(k),sbound*sigma(k));
            q.strength{k} = q.strength{k}/sum(q.strength{k}) *nactual(k);
        end
    end
% normalize strength to prescribed Neff
    ratio = (prescribedNeff/sum(prescribedNeff)) ./ (nactual/sum(nactual));
    for k=1:subregion
        q.strength{k} = s * q.strength{k} * ratio(k);
    end
    if sum(q.isub) == 0, return; end
    % recalibrate theta
    q.theta = prescribedTheta;
    %if subregion == 1
    %    q.theta = prescribedTheta;
    %end
    %if subregion == 2
    %    %if norm(peak(1,:)-peak(2,:)) > mean(sigma)*sqrt(2*log(2));
    %        q.theta = atan(diff(peak(:,2))/diff(peak(:,1)));
    %        if isnan(q.theta)
    %            q.theta = pi/2;
    %        end
    %    %else
    %        q.theta = prescribedTheta;
    %    %end
    %end
    %if subregion > 2
    %    f = fit(peak(:,1), peak(:,2),'poly1');
    %    q.theta = atan(f.p1); 
    %end
    %% assign on off
    if subregion == 1
        if rand < pOn
            q.nLGNon = q.isub(1);
        else
            q.v1Map{1} = -q.v1Map{1};
            q.nLGNoff = q.isub(1);
            q.isub(1) = -q.isub(1);
        end
    else
        if rand < 0.5
            assign = 1;
        else
            assign = -1;
        end

        for k=1:subregion
            if assign > 0
                q.nLGNon = q.nLGNon + q.isub(k);
                assign = -1;
            else
                q.v1Map{k} = -q.v1Map{k};
                q.nLGNoff = q.nLGNoff + q.isub(k);
                q.isub(k) = -q.isub(k);
                assign = 1;
            end
        end       
    end
    %% output
    q.sp = 0;
    for k=1:subregion
        if k > 1 && abs(q.isub(k))>0 && abs(q.isub(k-1))>0
            q.sp = q.sp + 2*norm(peak(k,:) - peak(k-1,:));
        end
    end
    if subregion > 1
        q.sp = q.sp/(subregion-1);
    else
        q.sp =  2*2*sqrt(2*log(2))*sigma;
    end
    if isnan(q.sp)
        error('sp nan');
    end
    if example~=0
        f = figure;
        mStrength = zeros(1,subregion);
        hold on
        ymax = -inf;
        ymin = inf;
        xmax = ymax;
        xmin = ymin;
        for k = 1:subregion       
            xxx2 = rfx(peak(k,1),sigmb(k),sigma(k),gtheta(k),ttt,FWHM,FWHM)*180/pi;         
            yyy2 = rfy(peak(k,2),sigmb(k),sigma(k),gtheta(k),ttt,FWHM,FWHM)*180/pi;
            xxx1 = rfx(peak(k,1),sigmb(k),sigma(k),gtheta(k),ttt,HWHM,HWHM)*180/pi;         
            yyy1 = rfy(peak(k,2),sigmb(k),sigma(k),gtheta(k),ttt,HWHM,HWHM)*180/pi;
            xxx0 = rfx(peak(k,1),sigmb(k),sigma(k),gtheta(k),ttt,boundb(k),bounda(k))*180/pi;         
            yyy0 = rfy(peak(k,2),sigmb(k),sigma(k),gtheta(k),ttt,boundb(k),bounda(k))*180/pi;
            ymin = min(min(yyy2),ymin);
            ymax = max(max(yyy2),ymax);
            xmin = min(min(xxx2),xmin);
            xmax = max(max(xxx2),xmax);
            if sum(q.v1Map{k}) > 0
                plot(xxx2,yyy2,'--r','LineWidth',1);
                plot(xxx1,yyy1,':r','LineWidth',1); 
                plot(xxx0,yyy0,'-r','LineWidth',1); 
%                 plot3(dLGNpos(q.v1Map{k},1), dLGNpos(q.v1Map{k},2),q.strength{k}*5,'*r','MarkerSize',4.0);
                plot(dLGNpos(q.v1Map{k},1), dLGNpos(q.v1Map{k},2),'*r','MarkerSize',4.0);
                plot(peak(k,1)*180/pi,peak(k,2)*180/pi,'sr');
                plot(opeak(k,1)*180/pi,opeak(k,2)*180/pi,'^r');
                x = [min(xxx2),max(xxx2)];
                plot(x,tan(gtheta(k))*(x-peak(k,1)*180/pi)+peak(k,2)*180/pi,'--r');
                
            else
                plot(xxx2,yyy2,'--b','LineWidth',1);
                plot(xxx1,yyy1,':b','LineWidth',1); 
                plot(xxx0,yyy0,'-b','LineWidth',1); 
%                 plot3(dLGNpos(-q.v1Map{k},1), dLGNpos(-q.v1Map{k},2),q.strength{k}*5,'ob','MarkerSize',4.0);
                plot(dLGNpos(-q.v1Map{k},1), dLGNpos(-q.v1Map{k},2),'ob','MarkerSize',4.0);
                plot(peak(k,1)*180/pi,peak(k,2)*180/pi,'sb');
                plot(opeak(k,1)*180/pi,opeak(k,2)*180/pi,'^b');
                x = [min(xxx2),max(xxx2)];
                plot(x,tan(gtheta(k))*(x-peak(k,1)*180/pi)+peak(k,2)*180/pi,'--b');
            end
            mStrength(k) = sum(q.strength{k});
        end
        ylim([ymin,ymax]);
        xlim([xmin,xmax]);
        dxy = rlgn*180/pi;
        edge = 1;
        xrange = xlim; x = (xrange(1)-edge*dxy):(dxy/24):(xrange(2)+edge*dxy);
        yrange = ylim; y = (yrange(1)-edge*dxy):(dxy/24):(yrange(2)+edge*dxy);   

%         zmax = 0;
        for k = 1:subregion
            if sum(q.v1Map{k}) > 0
                [X, Y, Z] = multiLGNspatialKernel(dLGNpos(q.v1Map{k},1), dLGNpos(q.v1Map{k},2),x,y,q.isub(k),q.strength{k},200);
                contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','m','LineWidth',2);
                newZmax = max(max(Z));
                contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','m','LineWidth',2);
%                 mesh(X,Y,Z,'FaceColor','None','EdgeColor','r');
            else
                [X, Y, Z] = multiLGNspatialKernel(dLGNpos(-q.v1Map{k},1), dLGNpos(-q.v1Map{k},2),x,y,-q.isub(k),q.strength{k},200);
                contour(X,Y,Z,[0,0],'LineStyle','--','LineColor','g','LineWidth',2);
                newZmax = max(max(Z));
                contour(X,Y,Z,[newZmax/2,newZmax/2],'LineStyle',':','LineColor','g','LineWidth',2);
%                 mesh(X,Y,Z,'FaceColor','None','EdgeColor','b');
            end
            
%             if  newZmax > zmax
%                 zmax = newZmax;
%             end
        end
        xlim([min(x),max(x)]);
        ylim([min(y),max(y)]);
%         zlim([0,inf]);
        grid on
        axis equal
        if q.theta < 0 
            q.theta = q.theta + pi;
        end
        if q.theta >= pi
            q.theta = q.theta - pi;
        end
		if prescribedTheta >= pi/2
			oPrescribedTheta = prescribedTheta - pi/2;
		else
			oPrescribedTheta = prescribedTheta + pi/2;
		end
		if q.theta >= pi/2
			otheta = q.theta - pi/2;
		else
			otheta = q.theta + pi/2;
		end
        title({['drifting/grating angle: ',num2str(prescribedTheta*180/pi,'%3.1f'),'^o, ',num2str(oPrescribedTheta*180/pi,'%3.1f'),'^o'],...
			['recalibrated: ',num2str(q.theta*180/pi,'%3.1f'),'^o, ',num2str(otheta*180/pi,'%3.1f'),'^o'],...
            ['prescribed Neff: ',num2str(round(prescribedNeff))],...
            ['get: ',num2str(q.isub')],...
            ['mStrength: ', num2str(mStrength,'% 2.1f')],...
            type});
                %    'strength are normalized to prescribed portion',...
        xlabel('degree');%, LGN rf radius:', num2str(rlgn*180/pi)]);
        ylabel('degree');

%         campos([xrange(1),yrange(1),zmax*1.5]);
        if ~isempty(format)
            set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(format,'fig')
                savefig(f,['lgn2v1map/',suffix,'/',type,'.',format]);
            else
                print(f,['lgn2v1map/',suffix,'/',type,'.',format],printDriver,dpi);
            end
        end
    end
    q.nLGN = q.nLGNon + q.nLGNoff;
    q.center = mean(LGNpos(q.pick,:),1);
    if subregion == 2 
        if norm(peak(1,:)-peak(2,:)) < mean(sigma)*sqrt(2*log(2))*Yidx*2;
            q.iORF = true;
        end
    end
end
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
%     for i = 1:n
%         Z = Z + 1*(Aa/(siga^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/siga^2/2)...
%             - Ab/(sigb^2*2*pi) * exp(-((X-posx(i)).^2+(Y-posy(i)).^2)/sigb^2/2));
%     end
    Z = Z/scale;
end
