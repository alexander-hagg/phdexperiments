%% Add to path
%workdir = getenv('PBS_O_WORKDIR');
if ~exist('workdir','var') || isempty(workdir);workdir='/home/alex/experiments';end
addpath(genpath([workdir '/utils']));
addpath(genpath([workdir '/clustering']));
addpath(genpath([workdir '/modules']));
rmpath(genpath([workdir '/modules/ffdSail/domains']));
addpath(genpath([workdir '/metrics']));
addpath(genpath([workdir '/tests']));
addpath(genpath([workdir '/dimensionalityreduction']));

if strcmp(domainname,'PARSEC'); addpath(genpath([workdir '/modules/ffdSail/domains/parsec'])); end
if strcmp(domainname,'FOILFFD'); addpath(genpath([workdir '/modules/ffdSail/domains/foilFFD'])); end

%% Parallelization
myCluster = parcluster('local');
numWorkers = str2num(getenv('PBS_NUM_PPN'));
if isempty(numWorkers);numWorkers = 4;end % is local    
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

