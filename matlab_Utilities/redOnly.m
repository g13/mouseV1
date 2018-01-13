function c=redOnly(m)
if nargin<1,m = size(get(gcf,'colormap'),1);end
    r0 = 0.0;
% From [constant 1 1] to [constant 0 0]
    r = ones(m,1);
    b = linspace(r0,1,m);
    b = fliplr(b);
    b = b';
    g = b;
    c = [r g b];
end
