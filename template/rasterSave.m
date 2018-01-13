function rasterSave(theme,lgnfile,eps,itheta,processed)
if nargin < 5
    processed = true;
end
load(lgnfile,'p');
epsStr = num2str(eps,'%04d');
thetaStr = num2str(itheta,'%02d'); 
 
if processed
    DIR = [theme,'/spike_wise/',epsStr,'/']; 
    [tspI,l] = readSpikes(DIR,p.nv1,[thetaStr,'-spikes.dat']);
else 
    DIR = [theme,'/',epsStr,'/',thetaStr,'/']; 
    [tspI,l] = readSpikes(DIR,p.nv1,['spikes.dat']);
end
save([theme,'-',epsStr,'-',thetaStr,'-rasterStruct'],'tspI','l');
end
