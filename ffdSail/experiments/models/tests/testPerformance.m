clear all; clc;

systemInit;
d = ffd_Domain;
p = sail;
numSamples = [50, 250, 500, 1000, 2000];
[modelPars,modelNames] = models_LoadCfg(d);
numReps = 1;
datafile = 'samples_50-2000.mat';
outpath = '/scratch/ahagg2s/performance/';
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
for rep=1:1%numReps
    outStr = ['Replicate ' int2str(rep) '\n'];
    fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
    for sampleSet=1:size(observation,1)
        outStr = [char(9) 'Sampleset ' int2str(sampleSet) '\n'];
        fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
        cat{sampleSet}(rep,:,:) = feval(d.categorize, testObservation{sampleSet,rep}, d);
        for modelIndex=1:length(modelPars)
            outStr = [char(9) char(9) 'Model: ' modelNames{modelIndex}];
            fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr);fclose(fid);
            cfg = feval(['params' modelPars{modelIndex}{1}] ,modelPars{modelIndex}{2:end});
            tic;
            models{sampleSet,modelIndex,rep} = feval(['train' modelPars{modelIndex}{1}], observation{sampleSet,rep}, value{sampleSet,rep}(:,1), cfg);
            perfTrain(modelIndex,sampleSet) = toc;
            tic;
            pred{sampleSet}(modelIndex,rep,:,:) = feval(['predict' modelPars{modelIndex}{1}], models{sampleSet,modelIndex,rep}, testObservation{sampleSet,rep});
            perfPredict(modelIndex,sampleSet) = toc;
            
            trues{sampleSet}(modelIndex,rep,:) = testValue{sampleSet,rep}(:,1)';
            errs{sampleSet}(modelIndex,rep,:) = pred{sampleSet}(modelIndex,rep,:,1)-trues{sampleSet}(modelIndex,rep,:);
            
            modelErrors = reshape(errs{sampleSet}(modelIndex,:,:),1,size(errs{sampleSet},2)*size(errs{sampleSet},3)*size(errs{sampleSet},4));
            truedata = reshape(trues{sampleSet}(modelIndex,:),1,size(trues{sampleSet},2)*size(errs{sampleSet},3));
            out1 = 100*prctile(modelErrors(:)/median(truedata(:)),25);
            out2 = 100*median(modelErrors(:))/median(truedata(:));
            out3 = 100*prctile(modelErrors(:)/median(truedata(:)),75);
            outStr = [char(9) '25Q - Med - 75Q error:    %.1f | %.1f | %.1f \n'];
            fid = fopen([outpath 'progress'],'a');fprintf(fid, outStr, out1, out2, out3);fclose(fid);
            
        end
    end
end

%% VISUALIZATION
clrs = parula(length(modelPars)+1);


% Get ranking scores for all samples
close all;
fig(1) = figure(1);hold off;
fig(2) = figure(2);hold off;

for modelIndex=1:length(modelPars)
    figure(1);
    semilogy(perfTrain(modelIndex,:)', 'x-','LineWidth',4, 'Color', clrs(modelIndex,:));
    hold on;
    figure(2);
    semilogy(perfPredict(modelIndex,:)', 'x-','LineWidth',4, 'Color', clrs(modelIndex,:));
    hold on;
end

for i=1:2
    figure(i);
    legend(string(modelNames),'Location','SouthEast');
    ax = gca;
    ax.XTick = 1:length(numSamples);
    if i==1;ax.XTickLabels = string(numSamples);end
    if i==2;ax.XTickLabels = string(numSamples/5);end
    %axis([0.9 0.1+size(observation,1) 50 100]);
    xlabel('Number of Samples');
    ylabel('Time taken (s)');
    if i==1;title('Training Performance');end
    if i==2;title('Prediction Performance');end
    grid on;
end

%%
save_figures(fig, outpath, ['performance:'], 24, [14 8]);

%%
save([outpath 'results.mat']);