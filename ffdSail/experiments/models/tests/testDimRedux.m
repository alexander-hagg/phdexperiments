clear all; clc;

systemInit;
d = ffd_Domain;
p = sail;

numSamples = [50, 250, 500, 1000];
cfg = paramsBootstrap(10);
pcaMin = [NaN 0.005 0.01 0.02];
numReps = 5;

disp('Initializing Sample Sets');
for ss=1:length(numSamples)
    disp(['Set ' int2str(ss) '/' int2str(length(numSamples)) ': ' int2str(numSamples(ss)) ' samples']);
    disp(['Train Set']);
    p.nInitialSamples = numSamples(ss);
    [observation{ss}, value{ss}] = initialSampling(d,p.nInitialSamples);
    
    % Normalization
    [observation{ss},psx{ss}] = mapminmax(observation{ss}');
    observation{ss} = observation{ss}';
    [value{ss},pst{ss}] = mapminmax(value{ss}');
    value{ss} = value{ss}';
    
    disp(['Test Set (20% of Training Set Size)']);
    p.nInitialSamples = numSamples(ss)/5;
    [testObservation{ss}, testValue{ss}] = initialSampling(d,p.nInitialSamples);
    testObservation{ss} = mapminmax('apply',testObservation{ss}',psx{ss});
    testObservation{ss} = testObservation{ss}';
    testValue{ss} = mapminmax('apply',testValue{ss}',pst{ss});
    testValue{ss} = testValue{ss}';
end

%%
clear pred trues errs models

for replicate=1:numReps
    disp(['Replicate ' int2str(replicate)]);
    for sampleSet=1:length(observation)
        disp(['Sampleset ' int2str(sampleSet)]);
        for pcaInd=1:length(pcaMin)
            disp(['pca min variance: ' num2str(pcaMin(pcaInd))]);
            modelcfg.cfg = paramsHSM(size(observation{sampleSet},1),cfg);
            modelcfg.cfg.experts.pca_min_variance = pcaMin(pcaInd);
            models{sampleSet,pcaInd,replicate} = trainHSM(observation{sampleSet}', value{sampleSet}(:,1), modelcfg);
            
            % Get dimensionality
            if pcaInd==1
                thesedims = [10 10 10 10 10];
            else
                nodes = models{sampleSet,pcaInd,replicate}.hierarchy.Node;
                for node=1:length(nodes)
                    thesedims(node) = nodes{node}{1}.cfg.pca.yrows;
                end
            end
            dims{sampleSet}(pcaInd,replicate,:) = thesedims;
            
            pred{sampleSet}(pcaInd,replicate,:) = predictHSM(testObservation{sampleSet}',models{sampleSet,pcaInd,replicate});
            trues{sampleSet}(pcaInd,replicate,:) = testValue{sampleSet}(:,1)';
            errs{sampleSet}(pcaInd,replicate,:) = pred{sampleSet}(pcaInd,replicate,:)-trues{sampleSet}(pcaInd,replicate,:);
            
            fig(sampleSet) = figure(sampleSet);
            subplot(2,1,1);
            data = squeeze(reshape(errs{sampleSet},size(errs{sampleSet},1),size(errs{sampleSet},2)*size(errs{sampleSet},3)))';
            boxplot(data);
            axis([0.5 length(pcaMin)+0.5 -1 1]);
            grid on;ax = gca;
            labels = string(pcaMin);labels(1) = 'NaN';ax.XTickLabel = labels;
            title([int2str(numSamples(sampleSet)) ' Samples']);
            ylabel('Error');xlabel('Max. Dimension Variance');
            drawnow;
            subplot(2,1,2);
            data = squeeze(reshape(dims{sampleSet},size(dims{sampleSet},1),size(dims{sampleSet},2)*size(dims{sampleSet},3)))';
            boxplot(data);
            axis([0.5 length(pcaMin)+0.5 0 10.5]);
            grid on;ax = gca;
            labels = string(pcaMin);labels(1) = 'NaN';ax.XTickLabel = labels;
            title([int2str(numSamples(sampleSet)) ' Samples']);
            ylabel('Num Dimensions Used');xlabel('Max. Dimension Variance');
            drawnow;
            
        end
    end
end

%%

% Get ranking scores
for pcaInd=1:length(pcaMin)
    for sampleSet=1:length(observation)
        
        pr{sampleSet}(:,:) = squeeze(pred{sampleSet}(pcaInd,:,:));
        cmpmat = repmat(pr{sampleSet}(:),1,length(pr{sampleSet}(:)));
        PRANK = cmpmat<pr{sampleSet}(:)';
        tr{sampleSet}(:,:) = squeeze(trues{sampleSet}(pcaInd,:,:));
        cmpmat = repmat(tr{sampleSet}(:),1,length(tr{sampleSet}(:)));
        TRANK = cmpmat<pr{sampleSet}(:)';
        validRankings = sum(PRANK(:)&TRANK(:)) + sum(~PRANK(:)&~TRANK(:));
        relValidRankings(pcaInd,sampleSet) = 100*validRankings'/(size(PRANK,1)*size(PRANK,2));
    end
end

fig(5) = figure(5);
plot(relValidRankings','LineWidth',2);
legend(labels,'Location','SouthEast');
ax = gca;
ax.XTick = 1:length(numSamples);
ax.XTickLabels = string(numSamples);
xlabel('Number of Samples');
ylabel('Ranking Accuracy (%)');
grid on;

%%
save_figures(fig, './', ['dimredux_error_'], 24, [12 8]);