% Get domain
d = mirror_Domain('nCases',1);
% Get shape
parameterLimit = 0.15;
mutation = parameterLimit*randn(1,45);
mutation(mutation<-parameterLimit) = -parameterLimit;
mutation(mutation>parameterLimit) = parameterLimit;
observation = 0.5+mutation;


%% Run OpenFoam
d.openFoamFolder = ['/scratch/ahagg2s/sailCFD/'];
d.openFoamTemplate = [d.runFolder '/ffdSail/domains/' d.namebase '/pe/ofTemplates/local_short/'];
%for iCase = 1:d.nCases
%    system(['rm -rf ' d.openFoamFolder 'case' int2str(iCase)]);
%    system(['mkdir ' d.openFoamFolder 'case' int2str(iCase)]);
%    system(['cp -r ' d.openFoamTemplate '* ' d.openFoamFolder 'case' int2str(iCase)]); 
%end

%system(['cd ' d.openFoamFolder ]);
%system(['./startCaseRunners.sh']);


%%
[value,timings] = mirror_PreciseEvaluate(observation, d);

disp(timings)
