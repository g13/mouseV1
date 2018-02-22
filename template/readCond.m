function [ipA,ioA,gtot,ge,gesc,gi,glgn,glgnsc,vs,vssc] = readCond(DIR,ntheta,n)
    filepath = [DIR, '/', 'cv.dat'];
    fid = fopen(filepath);
    Data = fread(fid, [n, 16],'double');
    ipA = nearest(Data(:,15));
    ioA = nearest(Data(:,16));
    fclose(fid);

    filepath = [DIR, '/', 'intra.dat'];
    fid = fopen(filepath);
    gtot =   fread(fid, [n, 2*ntheta+1],'double');
             fread(fid, [n, 2*ntheta+1],'double');
    vs   =   fread(fid, [n, 2*ntheta+1],'double');
    glgn =   fread(fid, [n, 2*ntheta+1],'double');
    ge =     fread(fid, [n, 2*ntheta+1],'double');
    gi =     fread(fid, [n, 2*ntheta+1],'double');
             fread(fid, [n, 2*ntheta+1],'double');
             fread(fid, [n, 2*ntheta+1],'double');
             fread(fid, [n, 2*ntheta+1],'double');
             fread(fid, [n, 2*ntheta+1],'double');
    glgnsc = fread(fid, [n, 2*ntheta+1],'double');
             fread(fid, [n, 2*ntheta+1],'double');
    gesc   = fread(fid, [n, 2*ntheta+1],'double');
             fread(fid, [n, 2*ntheta+1],'double');
    vssc   = fread(fid, [n, 2*ntheta+1],'double');
    fclose(fid);
end
