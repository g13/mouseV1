function [m,h] = subMat(src,tar,p,logP)
    
    if length(src.neff) == 1
        neff = ones(tar.n,1) * src.neff;
    else
        neff = src.neff;
    end
    if nargin == 4
        fid = fopen('logNormalProfile.out','w+');
        if ~isempty(p.profile)
            if logP.spread
                minDis = min(logP.normDistance);
                refDis = p.YoudensIndex;
                pow = log(0.1)/log((exp(0.2)-1)/(exp(refDis)-1));
                ratio = logP.lR + ((exp(logP.normDistance-minDis)-1)./(exp(refDis)-1)).^pow;
                ratio(logP.normDistance > refDis) = 1;
                profiles = zeros(logP.nbins+2,tar.n);
            else
                [~,slist] = p.profile(neff(1),logP.nbins,logP.mu,logP.sigma,logP.mu,logP.nsig,logP.p2s,false);
                fprintf(fid,'%G \n',[1,slist/logP.p2s,logP.mu0/logP.p2s]');
                profiles =[1,slist,logP.mu];
            end
        else
            profiles=[1];
        end
    end
    m = zeros(src.n,tar.n,'int8');
    if isempty(find(src.neff > 0,1))
        return
    end
    if ~isfield(p,'rep')
        p.rep = 1;
    end
    if std(neff) == 0
        hroll = rand(neff(1),1);
    end
    if p.self
        srca = src.a;
    else
        srca = src.a*src.n/(src.n-1);
    end
    %HWHM2 =  (src.raxn^2 + tar.rden^2);
    HWHM2 =  (src.raxn + tar.rden)^2;
    hlx = src.lx/2;
    hly = src.ly/2;

    if tar.bound > 0 % cutoff boundary
        r2 = tar.bound^2*HWHM2;
    else
        r2 = hlx^2 + hly^2;
    end
    needCoMat = false;
    switch p.case
                % distance, RFcoef, Ori ...
        case 0  % uniform
            A = r2*pi;
            g = @(pp) ones(pp.npick,1);%neff*srca/A;
            specific = false;
        case 1  % gauss
            sig2 = tar.sig * HWHM2/(2*log(2));
            g = @(pp) exp(-0.5*pp.x2/sig2);
            specific = false; 
        case 2  % uniform,  gauss
            p.specificMat2 = p.specificMat.^2; 
            sigCoeff2 = tar.sigCoeff.^2;
        	g = @(pp) exp(-0.5*pp.c2./sigCoeff2);
            specific = true;
            needCoMat = true;
        case 3  % gauss,    gauss
            p.specificMat2 = p.specificMat.^2;
            sig2 = tar.sig * HWHM2/(2*log(2));
            sigCoeff2 = tar.sigCoeff.^2;
            g = @(pp) exp(-0.5*(pp.x2./sig2+pp.c2./sigCoeff2));
            specific = true;
            needCoMat = true;
        case 4  % uniform, uniform
            %p.specificMat = p.specificMat
            g = @(pp) ones(pp.npick,1);%neff*srca/A;
            specific = true;
        case 5 % uniform, gauss
            sig2 = HWHM2/(2*log(2));
            g = @(pp) exp(-0.5*pp.x2/sig2);
            specific = true; 
        otherwise
            'non-defined connecting strategy';
            return
    end
           
    % (i,j) => i innervate j
    m_logical = false(src.n,tar.n);

    dx = zeros(src.n,1);
    dy = zeros(src.n,1);
    id0 = 1:src.n;
    rep = neff.*p.rep;
    roll0 = rand(sum(rep),1);
    
    n = 0;
    for i=1:tar.n
        %x
        pick = abs(src.x-tar.x(i)) > hlx;
        dx(pick) = 2*mod(tar.x(i)+hlx,src.lx) - tar.x(i) - src.x(pick);
        %assert(sum(dx(pick) >= hlx + 1e-5)==0);
        dx(~pick) = src.x(~pick) - tar.x(i);
        %y
        pick = abs(src.y-tar.y(i)) > hly;
        dy(pick) = 2*mod(tar.y(i)+hly,src.ly) - tar.y(i) - src.y(pick);
        %assert(sum(dy(pick) >= hly + 1e-5)==0);
        dy(~pick) = src.y(~pick) - tar.y(i);

        d2 = dx.^2 + dy.^2;
        % boundary cutoff
        if tar.bound > 0
            pick = d2 < r2;
        else
            pick = true(src.n,1);
        end
        % self connection
        if ~p.self
            pick(i) = false;
        end
        npick = sum(pick);
        pp.npick = npick;
        assert(npick > neff(i));
        d2 = d2(pick); 
        id = id0(pick);
        
        pp.x2 = d2;
        if needCoMat
            pp.c2 = p.specificMat2(pick,i);
        end

        quasi_prob = g(pp);
        bounds = cumsum(quasi_prob);
        cap = bounds(npick);

        
        roll = roll0(n+1:n+rep(i))*cap;
        %roll = roll*ones(1,npick);
        %sg = roll-ones(rep(i),1)*bounds';
        %idhist = sum(diff([false(rep(i),1),sg<0],1,2));

        idhist = zeros(npick,1);
        for j=1:rep(i)
            tempID = find((roll(j)-bounds)<0,1);
            idhist(tempID) = idhist(tempID) + 1;
        end
        pickedID = id(idhist > 0);
        npicked = length(pickedID);

        if p.exact
            if npicked >= neff(i)
                if npicked > neff(i)
                    shuffled = randperm(npick);
                    [~,ind] = sort(idhist(shuffled),'descend');
                    pickedID = id(shuffled(ind(1:neff(i))));
                end
                npicked = neff(i);
            else
                disp(['increase rep for exact presynaptic input, ', num2str(i)]);
                disp(['available ',num2str(npick),', ',num2str(npicked),'picked, need ',num2str(neff(i))]);
            end   
        end

        if ~isempty(p.profile)
            if logP.spread
                [~ ,ind] = sort(p.specificMat(pickedID,i)); 
                [m(pickedID(ind),i),slist] = p.profile(npicked,logP.nbins,logP.mu,logP.sigma*ratio(i),logP.mu0*(1-ratio(i))+logP.mu*ratio(i),logP.nsig,logP.p2s,false);
                fprintf(fid,'%G \n',[1,slist/logP.p2s,logP.mu0/logP.p2s]');
                profiles(:,i) =[1;slist';logP.mu0];
            else
                if specific
                    [~ ,ind] = sort(p.specificMat(pickedID,i)); 
                    %[m(pickedID(ind),i),~] = p.profile(npicked,rand(npicked,1));
                    [m(pickedID(ind),i),~] = p.profile(npicked,logP.nbins,logP.mu,logP.sigma,logP.mu,logP.nsig,logP.p2s,false);
                else 
                    [ss, ~] = p.profile(npicked,rand(npicked,1));
                    m(pickedID,i) = shuffle(ss);
                end
            end
        else
            m_logical(pickedID,i) = true;
        end
        n = n + rep(i);
    end

    if isempty(p.profile)
        m(m_logical) = 1;
    end
    if nargin == 4
        save('logNormalProfile.mat','profiles');
    end
    if nargin == 4
        fclose(fid);
        if logP.spread
            h = figure;
            subplot(2,1,1);
            ll = 51;
            edges = linspace(0,ceil(max(logP.normDistance)*100)/100,ll);
            [Counts,~,ind] = histcounts(logP.normDistance,edges);
            EPSPslice = zeros(ll-1,3);
            for i = 1:tar.n
                j = ind(i);
                neighbor = m(:,i)>0;
                tmp = profiles(m(neighbor,i),i);
                EPSPslice(j,1) = EPSPslice(j,1) + mean(tmp(:));
                EPSPslice(j,2) = EPSPslice(j,2) + min(tmp(:));
                EPSPslice(j,3) = EPSPslice(j,3) + max(tmp(:));
            end
            for j = 1:(ll-1)
                nind = sum(ind==j);
                EPSPslice(j,:) = EPSPslice(j,:)./nind;
            end
            xxx = (edges(1:ll-1) + edges(2:ll))./2;
            [ax,h1,h2] = plotyy(xxx,EPSPslice(:,1),xxx,Counts);
            h1.Color='k';
            h1.LineStyle='-';
            hold(ax(1));
            plot(ax(1),xxx,EPSPslice(:,2),':k');
            plot(ax(1),xxx,EPSPslice(:,3),':k');
            ylim(ax(1),[0,inf]);
            ylim(ax(2),[0,inf]);
            subplot(2,1,2);
            histogram(logP.normDistance,edges);
            xlabel('normalized Distance bwt On-Off Subregion');
        end
    end
end
