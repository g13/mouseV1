function [strengthID, slist] = logNormalProfile(n,nbins,mu,sigma,mu0,nsig,p2s,testProfile)
% setup1:ff logNormalProfile_beta(400,rand(400,1),30,0.28,0.5,true)
% Sen et al (2005) rats
% sig = 0.9355; % mean logrithmized
% m = -0.702; % std logrithmized
% Cossel et al (2015) mouse V1
% mu = 0.45; % mean non-logrithmized
% sigma = 0.68; % std non-logrithmized
% 
%   log(non-log)
%   m = log(mu/sqrt(1+sigma^2/mu^2);
%   sig = sqrt(log(1+sigma^2/mu^2));
%
%   non-log(log)
%   mu = exp(m + sig^2/2);
%   sigma = sqrt(exp(sigma^2)-1)*mu;

%             sigma = 0.68; % logrithmized std
%             mu = 0.45; % logrithmized mean
    if nargin < 7
        testProfile = false;
    end
    m = log(mu/sqrt(1+sigma^2/mu^2));
    sig = sqrt(log(1+sigma^2/mu^2));
%     nsig = max(ceil(sigma/mu),nsig);   
    strengthID = zeros(n,1);
    A = 2/erfc(-(log(nsig*sig)-m)/(sig*sqrt(2)));
    logNormalCDF = @(x,m,sig) 1/2*erfc(-(log(x)-m)/(sig*sqrt(2)))*A;
    
    %slist = zeros(nbins,1);
    %nbins0 = round(0.6*nbins);
    %x0 = linspace(m-nsig*sig,0,nbins0);
    %slist(1:nbins0) = exp(x0);
    %tmp = linspace(1,mu+nsig*sigma,nbins-nbins0+1);
    %slist(nbins0+1:nbins) = tmp(2:end);
    x0max = min(log(4.0),m+nsig*sig);
    x0min = max(log(5e-3),m-nsig*sig);
%     x0max = m + nsig*sig;
    x0 = exp(linspace(x0min,x0max,nbins));
    slist = x0;
%     x00 = [0,x0];
%     slist = (x00(1:nbins) + x00(2:nbins+1))/2;

    %slist = linspace(exp(m-nsig*sig),exp(m+nsig*sig),nbins);
        
%     x0max = 4.0;
%     x0min = 1e-2;
% %     x0max = m + nsig*sig;
%     x0 = linspace(x0min,x0max,nbins);
% %     x00 = [0,x0];
%     slist = x0;

    
    cuts = logNormalCDF(x0,m,sig);
    cuts = cuts/max(cuts);
%     disp(max(cuts));
    csnpick = round(cuts*n);
    npick = diff([0,csnpick]);
    last = 0;
    for i = 1:nbins
        if npick(i) > 0
            if last+npick(i) > n
               npick(i) = n - last; 
            end
            ipick = last+(1:npick(i));
            strengthID(ipick) = i;
            last = last + npick(i);
        end
        if last == n
            break;
        end
    end
    if testProfile
        %% testing
        figure;
        logNormal = @(x,m,sig) 1./(x*sqrt(2*pi)*sig).*exp(-(log(x)-m).^2/(2*sig^2));
        pick = npick>0;
        x = slist(pick);
        y = npick(pick);
        nx = sum(pick);
        dx = zeros(1,nx);
        for i = 1:nx
            switch i
                case 1
                    dx(i) = x(i+1)-x(i);
                case nx
                    dx(i) = x(i) - x(i-1);
                otherwise
                    dx(i) = 0.5*(x(i+1)-x(i-1));
            end
            
        end
        yprob = logNormal(x,m,sig).*dx;
        yprob = yprob./sum(yprob);

        subplot(2,1,1);
        semilogx(x,yprob*n);
        xlabel('PSP');
        ylabel('#');
        hold on

        semilogx(x,y);
        ylim([0,inf]);
        legend({'analytical','data'});
        meanx = sum(x.*yprob);
        stdx = sqrt(sum(x.^2.*yprob) - meanx^2);
        title({[num2str(meanx),'\pm',num2str(stdx)],['x0=',num2str(min(x0)),'~',num2str(max(x0))],...
                ['x=',num2str(min(x)),'~',num2str(max(x))]})

        subplot(2,1,2);
        plot(x,yprob*n);
        hold on

        plot(x,y);
        ylim([0,inf]);
        text(x(end),y(end)*2.0,num2str(y(end)));
        ss = slist(strengthID);
        means = mean(ss);
        stds = std(ss);
%         title({[num2str(means),'\pm',num2str(stds)],...
%             [num2str(mu),'\pm',num2str(sigma)]})
        
        yyprob = y./sum(y);
        assert(sum(y) == n);
        means0 = sum(x.*yyprob);
        assert(abs(means0-means)<1e-14);
        rt = mu0/means;
        xx = x*rt;
        plot(xx,y);
        edges = 0:0.2:5;
%         edges=linspace(slist(1),slist(end),nbins);
%         dedges=edges(2)-edges(1);
%         edges=[edges(1)-dedges/2, edges+dedges/2];
        histogram(ss*rt,edges);
        xlim([edges(1),inf]);
        meanss = sum(xx.*yyprob);
        stdss = sqrt(sum(xx.^2.*yyprob)-meanss^2);
        title({[num2str(means),'\pm',num2str(stds)],[num2str(meanss),'\pm',num2str(stdss)],...
            [num2str(mu),'\pm',num2str(sigma)],[num2str(nsig),'xnsig']})
        %
        slist = slist*rt;
        disp(['median=',num2str(median(ss))]);
    else
        if mu ~= mu0
            slist = (mu0/mean(slist(strengthID)))*slist; 
        end
    end
    strengthID = flipud(strengthID)+1;
end
