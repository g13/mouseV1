function c=phaseSkip(m)
if nargin<1,m = size(get(gcf,'colormap'),1);end
c0 = [0.8,1,0.8];
if m == 1
    c = c0;
else
    if (mod(m,2) == 0)
        % From [0 0 0] to [1 1 1], then [1 1 1] to [0 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        r = [r; flipud(r)];
        g = r;
        b = r;
    else
        % From [0 0 0] to [1 1 1] to [0 0 0];
        m1 = floor(m*0.5);
        r = (0:m1-1)'/max(m1,1);
        r = [r; 1; flipud(r)];
        g = r;
        b = r;
    end
    c1 = [r g b];
    c = [c0; c1];
end
end
