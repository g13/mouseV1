function [mTarSrc, neff0] = retroMat_beta(src,tar,p,reciprocal,mSrcTar)
    mTarSrc = zeros(src.n,tar.n,'int8');
    if isempty(find(src.neff > 0,1)) || reciprocal < 0
        return
    end
    if length(src.neff) == 1
        neff = ones(tar.n,1) * src.neff;
    else
        neff = src.neff;
    end   
    if ~isfield(p,'rep')
        p.rep = 1;
    end
    if ~isfield(p,'repRetro')
        p.repRetro = 1;
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
            %sig2 = HWHM2/(2*log(2));
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
            %sig2 = HWHM2/(2*log(2));
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
            %sig2 = HWHM2/(2*log(2));
            sig2 = tar.sig * HWHM2/(2*log(2));
            g = @(pp) exp(-0.5*pp.x2/sig2);
            specific = true; 
        otherwise
            'non-defined connecting strategy';
            return
    end
    switch p.rcase
                % distance, RFcoef, Ori ...
        case 0  % uniform
            A = r2*pi;
            %rg = @(pp) neff*srca/A;
            rg = @(pp) ones(pp.npick,1);
        case 1  % gauss
            %sig2 = HWHM2/(2*log(2));
            sig2 = tar.sig * HWHM2/(2*log(2));
            rg = @(pp) exp(-0.5*pp.x2/sig2);
        otherwise
            'non-defined connecting strategy';
            return
    end
    % (i,j) => i innervate j
    m_logical = false(src.n,tar.n);
    nRev = sum(mSrcTar>0,2);
    dx = zeros(src.n,1);
    dy = zeros(src.n,1);
    d2 = zeros(src.n,tar.n);

    neff0 = round(nRev * reciprocal);
%     neff0 = neff*reciprocal;
    rep = nRev*p.repRetro;
    % deal with retro connection
    pick0 = mSrcTar > 0;
    id0 = 1:src.n;
    roll0 = rand(sum(rep),1);
    n = 0;
    for i=1:tar.n
        % x
        pick = abs(src.x-tar.x(i)) > hlx;
        % enforce periodic boundary condition
        dx(pick) = 2*mod(tar.x(i)+hlx,src.lx) - tar.x(i) - src.x(pick);
        %assert(sum(dx(pick) >= hlx + 1e-5)==0);
        dx(~pick) = src.x(~pick) - tar.x(i);
        %y
        pick = abs(src.y-tar.y(i)) > hly;
        dy(pick) = 2*mod(tar.y(i)+hly,src.ly) - tar.y(i) - src.y(pick);
        %assert(sum(dy(pick) >= hly + 1e-5)==0);
        dy(~pick) = src.y(~pick) - tar.y(i);
        d2(:,i) = dx.^2 + dy.^2;
        % boundary cutoff
        if tar.bound > 0
            pick = d2(:,i) < r2 & pick0(i,:)';
        else
            pick = pick0(i,:);
        end
        npick = sum(pick);
        pp.npick = npick;
        id = id0(pick);
        pp.x2 = d2(pick,i);
        
        quasi_prob = rg(pp);
        bounds = cumsum(quasi_prob);
        cap = bounds(npick);

        roll = roll0(n+1:n+rep(i))*cap;
        %roll = roll*ones(1,npick);
        %sg = roll-ones(rep(i),1)*bounds';
        %idhist = sum(diff([false(rep(i),1),sg<0],1,2));

        idhist = zeros(nRev(i),1);
        for j=1:rep(i)
            tempID = find((roll(j)-bounds)<0,1);
            idhist(tempID) = idhist(tempID) + 1;
        end

        pickedID = id(idhist>0);
        npicked = length(pickedID);

        if p.exact
            if npicked >= neff0(i)
                if npicked > neff0(i)
                    shuffled = randperm(nRev(i));
                    [~, index] = sort(idhist(shuffled),'descend');
                    pickedID = id(shuffled(index(1:neff0(i))));
                end
%                     npicked = neff0(i);
            else
                disp('increase repRetro for exact reciprocal connection');
                disp(['available ',num2str(npick),', need ',num2str(neff(i))]);
            end
        end
        m_logical(pickedID,i) = true;
        n = n + rep(i);
    end

    pick0 = ~m_logical;
    neff0 = sum(m_logical)';
    
    neff = neff - neff0;
    if sum(neff<0) ~= 0
        disp(num2str(sum(neff<0)));
        error('neff < retroProb * Post.Neuron');
    end
   
    % deal with the remaining connection
    

    rep = neff.*p.rep;
    roll0 = rand(sum(rep),1);

    n = 0;
    for i=1:tar.n
        if neff(i) > 0
            % roll

            % boundary cutoff
            if tar.bound > 0
                pick = d2(:,i) < r2 & pick0(:,i);
            else
                pick = pick0(:,i);
            end
            npick = sum(pick);
            pp.npick = npick;
            assert(npick > neff(i));
            id = id0(pick);

            pp.x2 = d2(pick,i);
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
                if specific
                    [~ ,ind] = sort(p.specificMat(pickedID,i)); 
                    [m(pickedID(ind),i),~] = p.profile(npicked,rand(npicked,1));
                else 
                    [ss, ~] = p.profile(npicked,rand(npicked,1));
                    m(pickedID,i) = shuffle(ss);
                end
            else
                m_logical(pickedID,i) = true;
            end
            n = n + rep(i);
        end
    end
    if isempty(p.profile)
        mTarSrc(m_logical) = 1;
    end
end
