clear;
methods = {'UDHM','SEED'};


load('/scratch/ahagg2s/GECCO2019/planner_MUTATION/30UDHM.mat');
map(1) = runData{end}.acqMap;
clear runData

load('/scratch/ahagg2s/GECCO2019/planner_MUTATION/30SEED.mat');
map(2) = runData{3}.acqMap;

for methodID=1:length(map)    
    compareWith = [0 1 2 3];
    a = hist(map(methodID).ring1(:),compareWith);
    figure(methodID);
    viewMap(map(methodID).ring1,d);
    caxis([0 4]);colormap(parula(4));
    title(methods{methodID});
end

