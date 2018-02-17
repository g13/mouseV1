function cv = get_cv(d)
    assert(sum(d(:)<0)==0);
    n = size(d,1);
    theta = ((1:n)'./n*2*pi);
    repmatSize = size(d);
    repmatSize(1) = 1;

    theta = repmat(theta,repmatSize);
    value = sum(d.*exp(1j*theta));
    cv = squeeze(1-abs(value)./sum(d));
end
