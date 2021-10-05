%% Add to path
workdir = getenv('PBS_O_WORKDIR');
if ~exist('workdir','var') || isempty(workdir);workdir='/home/alex/ffdSail';end
addpath(genpath([workdir '/mapElites']));
addpath(genpath([workdir '/modules']));
addpath(genpath([workdir '/sail']));
addpath(genpath([workdir '/utils']));
addpath(genpath([workdir '/experiments']));
addpath(genpath([workdir '/visualization']));
addpath(genpath([workdir '/analysis']));


%addpath(genpath([workdir '/domains/foilFFD'])); 
if strcmp(domainname,'PARSEC'); addpath(genpath([workdir '/domains/parsec'])); end
if strcmp(domainname,'FOILFFD'); addpath(genpath([workdir '/domains/foilFFD'])); end
if strcmp(domainname,'MIRROR'); addpath(genpath([workdir '/domains/mirror'])); end

%% Parallelization
myCluster = parcluster('local');
% numWorkers = str2num(getenv('PBS_NUM_PPN'));
numWorkers = 4;
myCluster.NumWorkers = numWorkers;

% Create Temp Directory for Multithreading
tmpdir = getenv('TMPDIR');
if isempty(tmpdir);tmpdir='/tmp';end
myCluster.JobStorageLocation = tmpdir;
myCluster.HasSharedFilesystem = true;

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(myCluster, myCluster.NumWorkers);
end

%% Random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%% 
resultpath = getenv('EXPERIMENTPATH');
restart = getenv('RESTART');
if isempty(restart); restart=false;end
if ~exist('resultpath','var') || isempty(resultpath);resultpath='/scratch/ahagg2s/MISSING_EXPERIMENTPATH';end
mkdir(resultpath);

