clear;clc;
domainname = 'PARSEC';
cd('../..');systemInit;
resultpath = getenv('EXPERIMENTPATH');

%% Load results from mat files
resultFiles = {};
fsNodes = dir(resultpath);
fsNodes([fsNodes.isdir]) = [];

for f=1:length(fsNodes)
    resultFiles{end+1} = [resultpath '/' fsNodes(f).name];
end

%% Load result files
results={}; 

for iter=1:length(resultFiles)
   % Run maps through XFOIL maps
   disp(['Result file ' int2str(iter) '/' int2str(length(resultFiles))]);
   fname = resultFiles{iter};
   if strcmp(domainname,'PARSEC')
       dname = [workdir '/domains/parsec'];
   elseif strcmp(domainname,'FOILFFD')
       dname = [workdir '/domains/foilFFD'];
   end
   results{iter} = getGroundTruth(fname, dname, [], false);        
end

for f=1:length(fsNodes)
    output = results{f};
    p = output.p; d = output.d;
    save([resultpath '/' fsNodes(f).name], 'output', 'p', 'd', '-v7.3');
end

