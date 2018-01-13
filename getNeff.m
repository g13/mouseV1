%% Calculate mouse presynaptic neff
function [nEE, nEI, nIE, nII] = getNeff(e,i,layer)
den_V1 = 1.55e5; % per mm^3
if strcmp(layer,'2_3')
    depth_L2_3 = 0.19; % mm
    den = den_V1 * depth_L2_3;   % per mm^2, squeezed to plane.
end
if strcmp(layer,'4')
    depth_L4 = 0.11; % mm
    den = den_V1 * depth_L4;
end
% raxn and rden come in micro meter
% nEE = pi*(e.raxn^2 + e.rden^2)/1e6 *den* e.r*e.probe;
% nEI = pi*(i.raxn^2 + e.rden^2)/1e6 *den* i.r*e.probi;
% nIE = pi*(i.raxn^2 + e.rden^2)/1e6 *den* e.r*i.probe;
% nII = pi*(i.raxn^2 + i.rden^2)/1e6 *den* i.r*i.probi;

nEE = pi*(e.raxn + e.rden)^2/1e6 *den* e.r*e.probe;
nEI = pi*(i.raxn + e.rden)^2/1e6 *den* i.r*e.probi;
nIE = pi*(i.raxn + e.rden)^2/1e6 *den* e.r*i.probe;
nII = pi*(i.raxn + i.rden)^2/1e6 *den* i.r*i.probi;
end
