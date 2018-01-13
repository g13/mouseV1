
fdrlist = {'sample-5tf1-ie','sample-5t101-ie','sample-5t81-ie','sample-5t61-ie'}
theme = 'test50-ie';
pLabel = '\sigma of E<-I';
pTickLabel = {'flat','\sigma=1.0','\sigma=0.8','\sigma=0.6'};
lgnfile = '1xu-ndr305-50-s911.mat';
format = 'fig';
parameterTrace(fdrlist,pLabel,pTickLabel,lgnfile,format,theme)

fdrlist = {'sample-5tf5021','sample-5t5021','sample-5t8021','sample-5t6021'}
theme = 'test50';
pLabel = '\sigma of E<-I';
pTickLabel = {'flat','\sigma=1.0','\sigma=0.8','\sigma=0.6'};
lgnfile = '1xu-ndr305-50-s911.mat';
format = 'fig';
parameterTrace(fdrlist,pLabel,pTickLabel,lgnfile,format,theme)

fdrlist = {'sample-5tfe1','sample-5te1','sample-5t8e1','sample-5t6e1'}
theme = 'test50-nl';
pLabel = '\sigma of E<-I';
pTickLabel = {'flat','\sigma=1.0','\sigma=0.8','\sigma=0.6'};
lgnfile = '1xu-ndr305-50-s911.mat';
format = 'fig';
parameterTrace(fdrlist,pLabel,pTickLabel,lgnfile,format,theme)

fdrlist = {'sample-5tf50-rc','sample-5t102-rc','sample-5t8021','sample-5t60-rc'}
theme = 'test50-rc';
pLabel = '\sigma of E<-I';
pTickLabel = {'flat','\sigma=1.0','\sigma=0.8','\sigma=0.6'};
lgnfile = '1xu-ndr305-50-s911.mat';
format = 'fig';
parameterTrace(fdrlist,pLabel,pTickLabel,lgnfile,format,theme)

fdrlist = {'sample-5t8021'}
theme = 'test50-standard';
pLabel = '\sigma of E<-I';
pTickLabel = {'model'};
lgnfile = '1xu-ndr305-50-s911.mat';
format = 'fig';
parameterTrace(fdrlist,pLabel,pTickLabel,lgnfile,format,theme)
