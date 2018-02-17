function [spontRate,pkrate,rate] = readRate(DIR,ntheta,n,ndperiod)

filepath = [DIR, '/', 'frate.dat'];
fid = fopen(filepath);
rate = fread(fid, [n, 2*ntheta],'double');
fclose(fid);

filepath = [DIR, '/', 'cv.dat'];
fid = fopen(filepath);
Data = fread(fid, [n, 16],'double');
pkrate = Data(:,5);
spontRate = Data(:,9);
% neuronProperty.oorate = Data(:,11);
fclose(fid);

end
