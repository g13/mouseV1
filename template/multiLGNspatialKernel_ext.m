function [X,Y,Z] = multiLGNspatialKernel_ext(pos,x,y,n,s,p,scale)
    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));
    for i = 1:n
        Z = Z + s(i)*(p.Aa/(siga2*2*pi) * exp(-((X-pos(i,1)).^2+(Y-pos(i,2)).^2)/p.siga2/2)...
            - p.Ab/(p.sigb2*2*pi) * exp(-((X-pos(i,1)).^2+(Y-pos(i,2)).^2)/p.sigb2/2));
    end
    Z = Z/scale;
end
