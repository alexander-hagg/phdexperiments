%% Add to path
workdir = getenv('PBS_O_WORKDIR');
if ~exist('workdir','var') || isempty(workdir);workdir=pwd;end
addpath(genpath(workdir));
if strcmp(domainname,'PARSEC'); addpath(genpath([workdir '/domains/parsec'])); end
if strcmp(domainname,'FOILFFD'); addpath(genpath([workdir '/domains/foilFFD'])); end
if strcmp(domainname,'MIRROR'); addpath(genpath([workdir '/domains/mirror'])); end

%% Parallelization
myCluster = parcluster('local');
% numWorkers = str2num(getenv('PBS_NUM_PPN'));
numWorkers = 4;
myCluster.NumWorkers = 4;

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

