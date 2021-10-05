clear;clc;
numSelConcepts = 5;
domainname = 'FOILFFD';systemInit; xpfoldername = '/scratch/ahagg2s/acq_ffd_1'; shownRuns = [1];

load('data_prototypeselection_1.mat');

numReplicates = 10;

%%
for run=1:length(shownRuns)
    optima = reshape(data{1}.predMap(end/2).genes,625,10);
    optima(any(isnan(optima')),:) = [];
    
    labels = o{run}.estimatedLabels;
    
    
    for rep=1:numReplicates
        tree = fitctree(optima,labels); %view(tree,'mode','graph')
        predictThese = tree.predict(optima);
        prediction(run,rep,1,:) = predictThese;
        
        confMatInput1 = zeros(max(labels)+1,length(labels));
        for i=1:length(labels); confMatInput1(labels(i)+1,i) = 1;end
        confMatInput2 = zeros(max(labels)+1,length(labels));
        for i=1:length(labels); confMatInput2(predictThese(i)+1,i) = 1;end
        [c,~,~,per1] = confusion(confMatInput1,confMatInput2);
        
        confusionMatrix(run,rep,1,:) = sum(per1)/numel(unique(labels));
        
        treeEnsemble = TreeBagger(10,optima,labels, 'Method','classification');
        predictThese = treeEnsemble.predict(optima);
        predictThese = str2mat(predictThese);
        predictThese = str2num(predictThese);
        prediction(run,rep,2,:) = predictThese;
        
        confMatInput3 = zeros(max(labels)+1,length(labels));
        for i=1:length(labels); confMatInput3(predictThese(i)+1,i) = 1;end
        [c,~,~,per2] = confusion(confMatInput1,confMatInput3);
        
        confusionMatrix(run,rep,2,:) = sum(per2)/numel(unique(labels));
        
    end
    %view(treeEnsemble.Trees{2},'Mode','graph')
    
    
    
end

