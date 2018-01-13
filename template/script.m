theme =  'sample';
TYPE = 'M';
lgn = 'nd305-32';
eps = 1000;
ctheta = 6;
nthres = 20;
cthres = 3;
twindow = 5;
t1 = 21;
mdt = 5;
processed = true;
clusterFile = '';
%clusterOnly = false;
clusterOnly = true;

%clusterFile = ['Cluster-',theme,'-',num2str(nthres),'-',num2str(cthres),'.mat'];clusterOnly = false;

format = 'fig';
lgnfile = ['1xu-',lgn,'-s911.mat'];
conMatfile = [theme,'/conMat.mat'];
coMatfile = ['coMat-',lgn,'_',TYPE,'.mat'];
sfile = '';
Mov = false; %not implemented

clusterFile = rasterCluster(theme,lgnfile,eps,ctheta,conMatfile,coMatfile,sfile,clusterFile,nthres,cthres,Mov,format,processed,clusterOnly)
disp(clusterFile);

VeffCluster(theme,eps,ctheta,clusterFile,twindow,t1,format,mdt)
