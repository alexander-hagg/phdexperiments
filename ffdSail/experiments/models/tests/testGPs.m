clear all; clc;

systemInit;
d = ffd_Domain;
p = sail;
numSamples = [50, 250, 500, 1000];
% Get data
disp('Initializing Sample Sets');
datafile = 'samples_50-2000.mat';

[modelPars,modelNames] = models_LoadCfg(d);
numReps = 5;
disp(['Datafile ' datafile ' exists, loading...']);
load(datafile);
observation = {observation{3,:}};
value = {value{3,:}};
testObservation = {testObservation{3,:}};
testValue = {testValue{3,:}};
outpath = '/scratch/ahagg2s/custGP/';
mkdir(outpath);




%% Train and test models
for rep=1:numReps
    outStr = ['Replicate ' int2str(rep) '\n'];
    fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
    for sampleSet=1:size(observation,1)
        outStr = [char(9) 'Sampleset ' int2str(sampleSet) '\n'];
        fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
        catTrain{sampleSet}(rep,:,:) = feval(d.categorize, observation{sampleSet,rep}, d);
        cat{sampleSet}(rep,:,:) = feval(d.categorize, testObservation{sampleSet,rep}, d);
        for modelIndex=1:length(modelPars)
            outStr = [char(9) char(9) 'Model: ' modelNames{modelIndex} '\n'];
            fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
            cfg = feval(['params' modelPars{modelIndex}{1}] ,modelPars{modelIndex}{2:end});
            models{sampleSet,modelIndex,rep} = feval(['train' modelPars{modelIndex}{1}], observation{sampleSet,rep}, value{sampleSet,rep}(:,1), cfg);
            
            % Train error
            predTrain{sampleSet}(modelIndex,rep,:,:) = feval(['predict' modelPars{modelIndex}{1}], models{sampleSet,modelIndex,rep}, observation{sampleSet,rep});
            truesTrain{sampleSet}(modelIndex,rep,:) = value{sampleSet,rep}(:,1)';
            errsTrain{sampleSet}(modelIndex,rep,:) = predTrain{sampleSet}(modelIndex,rep,:,1)-truesTrain{sampleSet}(modelIndex,rep,:);
            
            modelErrors = reshape(errsTrain{sampleSet}(modelIndex,:,:),1,size(errsTrain{sampleSet},2)*size(errsTrain{sampleSet},3)*size(errsTrain{sampleSet},4));
            truedata = reshape(truesTrain{sampleSet}(modelIndex,:),1,size(truesTrain{sampleSet},2)*size(errsTrain{sampleSet},3));
            out1 = 100*median(modelErrors(:))/median(truedata(:));
            out2 = 100*abs(prctile(modelErrors(:)/median(truedata(:)),75)-prctile(modelErrors(:)/median(truedata(:)),25));
            outStr = [char(9) char(9) 'TRAIN MED - 25Q/75Q range perc. error:    %.1f | %.1f \n'];
            fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr, out1, out2);fclose(fid);
            
            % Test error            
            pred{sampleSet}(modelIndex,rep,:,:) = feval(['predict' modelPars{modelIndex}{1}], models{sampleSet,modelIndex,rep}, testObservation{sampleSet,rep});
            trues{sampleSet}(modelIndex,rep,:) = testValue{sampleSet,rep}(:,1)';
            errs{sampleSet}(modelIndex,rep,:) = pred{sampleSet}(modelIndex,rep,:,1)-trues{sampleSet}(modelIndex,rep,:);
            
            modelErrors = reshape(errs{sampleSet}(modelIndex,:,:),1,size(errs{sampleSet},2)*size(errs{sampleSet},3)*size(errs{sampleSet},4));
            truedata = reshape(trues{sampleSet}(modelIndex,:),1,size(trues{sampleSet},2)*size(errs{sampleSet},3));
            out1 = 100*median(modelErrors(:))/median(truedata(:));
            out2 = 100*abs(prctile(modelErrors(:)/median(truedata(:)),75)-prctile(modelErrors(:)/median(truedata(:)),25));
            outStr = [char(9) char(9) 'TEST  MED - 25Q/75Q range perc. error:    %.1f | %.1f \n'];
            fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr, out1, out2);fclose(fid);
        end
    end
end

%% VISUALIZATION
clrs = hsv(length(modelPars));
lstyle = {'-','-'};
lw = 4;

%%
figure(5);
subplot(2,2,1);hold off;
allTrues = squeeze(truesTrain{1}(1,:,:,1));
allTrues = allTrues(:);
[allTrues,id] = sort(allTrues);
plot(allTrues,'k-','LineWidth',2);
hold on;
t1 = squeeze(predTrain{1}(1,:,:,1));
t2 = squeeze(predTrain{1}(2,:,:,1));
t1 = t1(:);t2 = t2(:);
AllPreds = [t1(id)';t2(id)'];
plot(AllPreds','-','LineWidth',1);
axis([0 2500 -5 -2.5]);
legend(string({'Target',modelNames{:}}),'Location','NorthWest');

subplot(2,2,3);hold off;
boxplot((AllPreds-allTrues')');
grid on;

subplot(2,2,2);hold off;
allTrues = squeeze(trues{1}(1,:,:,1));
allTrues = allTrues(:);
[allTrues,id] = sort(allTrues);
plot(allTrues,'k-','LineWidth',2);
hold on;
t1 = squeeze(pred{1}(1,:,:,1));
t2 = squeeze(pred{1}(2,:,:,1));
t1 = t1(:);t2 = t2(:);
AllPreds = [t1(id)';t2(id)'];
plot(AllPreds','-','LineWidth',1);
axis([0 500 -5 -2.5]);
legend(string({'Target',modelNames{:}}),'Location','NorthWest');

subplot(2,2,4);hold off;
boxplot((AllPreds-allTrues')');
grid on;

%%
% 
% selection = 1:length(modelPars);
% % Get ranking scores for all samples
% fig(1) = figure(1);hold off;
% h{1} = viewRanking(predTrain, truesTrain, lstyle, clrs, lw, selection);
% xlabel('Number of Samples');ylabel('Ranking Error (%)');title('All comparisons (Train)');
% 
% % Get ranking scores for samples only within cells
% fig(2) = figure(2);hold off;
% for i=1:length(d.featureRes); edges{i} = linspace(0,1,d.featureRes(i)+1); end; %#ok<AGROW>
% h{2} = viewCellRanking(predTrain, truesTrain, catTrain, edges, lstyle, clrs, lw, selection);
% xlabel('Number of Samples');ylabel('Ranking Error (%)');title('Cell comparisons (Train)');
% 
% % Get ranking scores for all samples
% fig(3) = figure(3);hold off;
% h{3} = viewRanking(pred, trues, lstyle, clrs, lw, selection);
% xlabel('Number of Samples');ylabel('Ranking Error (%)');title('All comparisons (Test)');
% 
% % Get ranking scores for samples only within cells
% fig(4) = figure(4);hold off;
% for i=1:length(d.featureRes); edges{i} = linspace(0,1,d.featureRes(i)+1); end; %#ok<AGROW>
% h{4} = viewCellRanking(pred, trues, cat, edges, lstyle, clrs, lw, selection);
% xlabel('Number of Samples');ylabel('Ranking Error (%)');title('Cell comparisons (Test)');

% 
% for i=1:4
%     figure(i);
%     legend(h{i},string(modelNames),'Location','NorthWest');ax = gca;ax.XTick = 1:length(numSamples);ax.XTickLabels = string(numSamples/5);axis([0.9 0.1+size(observation,1) 0 50]);grid on;
% end

save_figures(fig, outpath, ['rank_error_'], 18, [14 8]);



%%
save([outpath 'results.mat'], '-v7.3');