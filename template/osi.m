function osi = osi(p,o)
    p(p<0) = 0;
    o(o<0) = 0;
    osi = (p-o)./(p+o);
    osi(isnan(osi)) = 0;
    osi(p+o==0) = 0;
end
