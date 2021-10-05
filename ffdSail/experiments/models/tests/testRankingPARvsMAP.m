clear all; clc;

systemInit;
d = ffd_Domain;
p = sail;
numSamples = [50, 250, 500, 1000, 2000];
[modelPars,modelNames] = models_LoadCfg(d);
numReps = 5;
datafile = 'samples_50-2000.mat';
outpath = '/scratch/ahagg2s/numSegments/';
mkdir(outpath);

% Get data
disp('Initializing Sample Sets');
if ~exist(datafile,'file')
    for rep=1:numReps
        disp(['Replicate: ' int2str(rep) '/' int2str(numReps)]);
        for ss=1:length(numSamples)
            disp(['Set ' int2str(ss) '/' int2str(length(numSamples)) ': ' int2str(numSamples(ss)) ' samples']);
            disp(['Train Set']);
            p.nInitialSamples = numSamples(ss);
            [observation{ss,rep}, value{ss,rep}] = initialSampling(d,p.nInitialSamples);
            
            disp(['Test Set (20% of Training Set Size)']);
            p.nInitialSamples = numSamples(ss)/5;
            [testObservation{ss,rep}, testValue{ss,rep}] = initialSampling(d,p.nInitialSamples);
        end
    end
    save(datafile, 'observation', 'value', 'testObservation', 'testValue');
else
    disp(['Datafile ' datafile ' exists, loading...']);
    load(datafile);
end


%% Train and test models
for rep=2:numReps
    outStr = ['Replicate ' int2str(rep) '\n'];
    fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
    for sampleSet=1:size(observation,1)
        outStr = [char(9) 'Sampleset ' int2str(sampleSet) '\n'];
        fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
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

%% Test models
for rep=1:size(pred{sampleSet},2)
    for sampleSet=1:size(observation,1)
        catTrain{sampleSet}(rep,:,:) = feval(d.categorize, observation{sampleSet,rep}, d);
        for modelIndex=1:length(modelPars)
            outStr = [char(9) char(9) 'Model: ' modelNames{modelIndex}];
            sprintf(outStr);
            predTrain{sampleSet}(modelIndex,rep,:,:) = feval(['predict' modelPars{modelIndex}{1}], models{sampleSet,modelIndex,rep}, observation{sampleSet,rep});
            truesTrain{sampleSet}(modelIndex,rep,:) = value{sampleSet,rep}(:,1)';
            errsTrain{sampleSet}(modelIndex,rep,:) = predTrain{sampleSet}(modelIndex,rep,:,1)-truesTrain{sampleSet}(modelIndex,rep,:);
            modelErrorsTrain = reshape(errsTrain{sampleSet}(modelIndex,:,:),1,size(errsTrain{sampleSet},2)*size(errsTrain{sampleSet},3)*size(errsTrain{sampleSet},4));
            truedataTrain = reshape(truesTrain{sampleSet}(modelIndex,:),1,size(truesTrain{sampleSet},2)*size(truesTrain{sampleSet},3));
        end
    end
end

%% VISUALIZATION
clrs = hsv(length(modelPars)+1);
lw = 2;

selection = [1,2,3,4,5;1,6,7,8,9];

for segWay=1:2
    
    % Get ranking scores for all samples
    fig((segWay-1)*4 + 1) = figure((segWay-1)*4 + 1);hold off;
    h{1} = viewRanking(predTrain, truesTrain, clrs(selection(segWay,:),:), lw, selection(segWay,:));
    xlabel('Number of Samples');ylabel('Ranking Error (%)');title('All comparisons (Train)');
    
    % Get ranking scores for samples only within cells
    fig((segWay-1)*4 + 2) = figure((segWay-1)*4 + 2);hold off;
    for i=1:length(d.featureRes); edges{i} = linspace(0,1,d.featureRes(i)+1); end; %#ok<AGROW>
    h{2} = viewCellRanking(predTrain, truesTrain, catTrain, edges, clrs(selection(segWay,:),:), lw, selection(segWay,:));
    xlabel('Number of Samples');ylabel('Ranking Error (%)');title('Cell comparisons (Train)');
    
    % Get ranking scores for all samples
    fig((segWay-1)*4 + 3) = figure((segWay-1)*4 + 3);hold off;
    h{3} = viewRanking(pred, trues, clrs(selection(segWay,:),:), lw, selection(segWay,:));
    xlabel('Number of Samples');ylabel('Ranking Error (%)');title('All comparisons (Test)');
    
    % Get ranking scores for samples only within cells
    fig((segWay-1)*4 + 4) = figure((segWay-1)*4 + 4);hold off;
    for i=1:length(d.featureRes); edges{i} = linspace(0,1,d.featureRes(i)+1); end; %#ok<AGROW>
    h{4} = viewCellRanking(pred, trues, cat, edges, clrs(selection(segWay,:),:), lw, selection(segWay,:));
    xlabel('Number of Samples');ylabel('Ranking Error (%)');title('Cell comparisons (Test)');
    
    for i=1:4
        figure((segWay-1)*4 + i);
        legend(h{i},string(modelNames(selection(segWay,:))),'Location','NorthWest');ax = gca;ax.XTick = 1:length(numSamples);ax.XTickLabels = string(numSamples/5);axis([0.9 0.1+size(observation,1) 0 25]);grid on;
    end
end

%%
save_figures(fig, outpath, ['rank_error_'], 18, [14 8]);

%%
save([outpath 'results.mat'], '-v7.3');