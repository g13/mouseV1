function [tick, n0] = autoAxis(x0,x1,n,lim,edge)
    common_ds = [2,2.5,5,10];
    if nargin <5
        edge = 1;
        if nargin<4
            lim = [-inf,inf];
            if nargin <3
                n = 6;
            end
        end
    end
    l = x1 - x0;
    d0 = l/n;
    e = ceil(log10(d0)-1);
    d = d0 * 10^(-e);
    di = abs(d - common_ds)./common_ds;
    [~,idx] = min(di);
    d = common_ds(idx) * 10^e;
    u0 = max(x0-edge*d, lim(1));
    u0 = floor(u0/d)*d + mod(u0,d); 
    u1 = min(x1+edge*d, lim(2));
    u1 = ceil(u1/d)*d + mod(u1,d); 
    tick = u0:d:u1;
    n0 = length(tick);
end
