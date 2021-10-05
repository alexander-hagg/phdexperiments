clear;clc;
load('data_prototypeselection_1.mat');
clear poolobj; %Remove old parallel pool object
domainname = 'FOILFFD';systemInit; xpfoldername = '/scratch/ahagg2s/acq_ffd_1'; 
shownRuns = [1];

for run=1:length(shownRuns)
    [experimentNames,data] = read_experiments(xpfoldername,'runsToShow', shownRuns(run));
    data{1}.d.fitnessPenalty_Area = true;
    data{1}.d.tmpdir = '/tmp/foil';
    tic;
    for c = 2:length(clusters)        
        o{run}.sail{c}.predMap = getGroundTruthMapArray(o{run}.sail{c}.predMap, data{1}.d, 'dummy', false);
    end
    toc
end





