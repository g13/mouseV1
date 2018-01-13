function c=blueOnly(m)
if nargin<1,m = size(get(gcf,'colormap'),1);end
    b0 = 0.0;
% From [1 1 constant] to [0 0 constant]
    b = ones(m,1);
    r = linspace(b0,1,m);
    r = fliplr(r);
    r = r';
    g = r;
    c = [r g b];
end
