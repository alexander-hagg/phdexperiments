clear all; clc;

systemInit;
d = ffd_Domain;
p = sail;
numSamples = [50, 250, 500, 1000, 2000];
[modelPars,modelNames] = models_LoadCfg(d);
numReps = 5;
datafile = 'samples_50-2000.mat';
outpath = '/scratch/ahagg2s/INvsINCAT/';
mkdir(outpath);

disp(['Datafile ' datafile ' exists, loading...']);
load(datafile);
observation = observation(1:4,:);value = value(1:4,:);testObservation = testObservation(1:4,:);testValue = testValue(1:4,:);

observation = [observation,observation,observation,observation,observation];
value = [value,value,value,value,value];
testObservation = [testObservation,testObservation,testObservation,testObservation,testObservation];
testValue = [testValue,testValue,testValue,testValue,testValue];

%% Train and test models
for rep=1:numReps
    outStr = ['Replicate ' int2str(rep) '\n'];
    fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
    for sampleSet=1:size(observation,1)
        outStr = [char(9) 'Sampleset ' int2str(sampleSet) '\n'];
        fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
        cat{sampleSet}(rep,:,:) = feval(d.categorize, testObservation{sampleSet,rep}, d);
        catTrain{sampleSet}(rep,:,:) = feval(d.categorize, observation{sampleSet,rep}, d);
        for modelIndex=1:length(modelPars)
            outStr = [char(9) char(9) 'Model: ' modelNames{modelIndex} '\n'];
            fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
            cfg = feval(['params' modelPars{modelIndex}{1}] ,modelPars{modelIndex}{2:end});
            
            if mod(modelIndex,2)==1
                input = observation{sampleSet,rep};
                inputTest = testObservation{sampleSet,rep};
            else
                input = [observation{sampleSet,rep} squeeze(catTrain{sampleSet}(rep,:,:))];
                inputTest = [testObservation{sampleSet,rep} squeeze(cat{sampleSet}(rep,:,:))];
            end
            
            models{sampleSet,modelIndex,rep} = feval(['train' modelPars{modelIndex}{1}], input, value{sampleSet,rep}(:,1), cfg);
            
            % Train error
            predTrain{sampleSet}(modelIndex,rep,:,:) = feval(['predict' modelPars{modelIndex}{1}], models{sampleSet,modelIndex,rep}, input);
            truesTrain{sampleSet}(modelIndex,rep,:) = value{sampleSet,rep}(:,1)';
            errsTrain{sampleSet}(modelIndex,rep,:) = predTrain{sampleSet}(modelIndex,rep,:,1)-truesTrain{sampleSet}(modelIndex,rep,:);
            
            modelErrors = reshape(errsTrain{sampleSet}(modelIndex,:,:),1,size(errsTrain{sampleSet},2)*size(errsTrain{sampleSet},3)*size(errsTrain{sampleSet},4));
            truedata = reshape(truesTrain{sampleSet}(modelIndex,:),1,size(truesTrain{sampleSet},2)*size(errsTrain{sampleSet},3));
            out1 = 100*median(modelErrors(:))/median(truedata(:));
            out2 = 100*abs(prctile(modelErrors(:)/median(truedata(:)),75)-prctile(modelErrors(:)/median(truedata(:)),25));
            outStr = [char(9) char(9) 'TRAIN MED - 25Q/75Q range perc. error:    %.1f | %.1f \n'];
            fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr, out1, out2);fclose(fid);
            
            % Test error
            pred{sampleSet}(modelIndex,rep,:,:) = feval(['predict' modelPars{modelIndex}{1}], models{sampleSet,modelIndex,rep}, inputTest);
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
clrs = repelem(hsv(length(modelPars)/2),2,1);
lstyle = {'-','--'};
lw = 4;
selection = 1:length(modelPars);
% Get ranking scores for all samples
fig(1) = figure(1);hold off;
h{1} = viewRanking(predTrain, truesTrain, lstyle, clrs, lw, selection);
xlabel('Number of Samples');ylabel('Ranking Error (%)');title('All comparisons (Train)');

% Get ranking scores for samples only within cells
fig(2) = figure(2);hold off;
for i=1:length(d.featureRes); edges{i} = linspace(0,1,d.featureRes(i)+1); end; %#ok<AGROW>
h{2} = viewCellRanking(predTrain, truesTrain, catTrain, edges, lstyle, clrs, lw, selection);
xlabel('Number of Samples');ylabel('Ranking Error (%)');title('Cell comparisons (Train)');

% Get ranking scores for all samples
fig(3) = figure(3);hold off;
h{3} = viewRanking(pred, trues, lstyle, clrs, lw, selection);
xlabel('Number of Samples');ylabel('Ranking Error (%)');title('All comparisons (Test)');

% Get ranking scores for samples only within cells
fig(4) = figure(4);hold off;
for i=1:length(d.featureRes); edges{i} = linspace(0,1,d.featureRes(i)+1); end; %#ok<AGROW>
h{4} = viewCellRanking(pred, trues, cat, edges, lstyle, clrs, lw, selection);
xlabel('Number of Samples');ylabel('Ranking Error (%)');title('Cell comparisons (Test)');

for i=1:4
    figure(i);
    legend(h{i},string(modelNames),'Location','NorthWest');ax = gca;ax.XTick = 1:length(numSamples);ax.XTickLabels = string(numSamples/5);axis([0.9 0.1+size(observation,1) 0 25]);grid on;
end


%
save_figures(fig, outpath, ['rank_error_'], 24, [14 8]);

%%
save([outpath 'results.mat'], '-v7.3');