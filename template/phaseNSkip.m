function c=phaseNSkip(m)
if nargin<1,m = size(get(gcf,'colormap'),1);end
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
    c = [r g b];
end
