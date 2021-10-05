% Hard Constraints (1)
% Similarity Mapping, (2)
% Can we learn categories 

% Categorization methods:
% - Non-linear: ANN, DT, RF
%
% 1.  Learn categorization
%       1a. Accuracy of categorization from original parameter space
%       1b. Accuracy of categorization from similarity space
% if 1b>1a
% 2.  Learn similarity space from original parameter space
%       2a. parameterized tSNE(X)
%       2b. GP(X->tSNE(X))

%% Load data
%clear; load('/scratch/ahagg2s/PPSN2018/SAILvsPRODUQD.mat','output');

inclZero = false;

trainRatio = 0.7;
valRatio = 0.15;
testRatio = 0.15;
for run=2:length(output)
    disp(['Run ' int2str(run)]);
    for iter=1:length(output{run}.data)
        disp(['Iter ' int2str(iter)]);
        % Get inputs and categories
        optima = output{run}.data{iter}.optima;
        latent = output{run}.data{iter}.latent;
        labels = output{run}.data{iter}.conceptLabels;
        if inclZero
            labels = labels + 1; % Rename zero label
        else
            noCat = labels==0; optima(noCat,:) = []; latent(noCat,:) = []; labels(noCat,:) = [];
        end
        
        Q = size(optima,1);
        allLabels = unique(labels);
        classlabels = (allLabels==labels')';
        % Train test split
        [trainInd,valInd,testInd] = dividerand(Q,trainRatio,valRatio,testRatio);
        
        %% Classification
        for inp = 1:2
            if inp==1
                inputs = optima;
                prototypes = output{run}.data{iter}.prototypes;
            elseif inp==2
                inputs = latent;
                prototypes = output{run}.data{iter}.prototypes;
                [~, index]=ismember(prototypes,optima,'rows');
                prototypes = latent(index,:);
            end
            clear outputs;
            % ANN
            % Create a Pattern Recognition Network
            hiddenLayerSize = 25;
            net = patternnet(hiddenLayerSize); net.trainParam.showWindow = 0;
            % Set up Division of Data for Training, Validation, Testing
            net.divideParam.trainRatio = trainRatio + valRatio;net.divideParam.valRatio = valRatio;net.divideParam.testRatio = 0;
            % Train the Network
            [net,tr] = train(net,inputs([trainInd,valInd],:)',classlabels([trainInd,valInd],:)');
            % Test the Network
            out = net(inputs(testInd,:)');
            [~,outputs(1,:)] = max(out);
            
            % DT
            Mdl = fitctree(inputs([trainInd,valInd],:),labels([trainInd,valInd]));
            outputs(2,:) = Mdl.predict(inputs(testInd,:));
            
            % RF
            iNumBags = 10;
            Mdl = TreeBagger(iNumBags,inputs([trainInd,valInd],:),labels([trainInd,valInd]),'OOBPred','On','Method','classification');
            out = predict(Mdl,inputs(testInd,:));
            outputs(3,:) = cellfun(@str2num, out(:,1:end))';
            
            % Distance to prototype
            distances = pdist2(inputs([testInd],:),prototypes);
            [~,outputs(4,:)] = min(distances');
            
            % Distance to class members
            compareto = [trainInd,valInd];
            distances = pdist2(inputs([testInd],:),inputs(compareto,:));
            [~,classMemberIDs] = min(distances');
            trainLabels = labels(compareto);
            outputs(5,:) = trainLabels(classMemberIDs);
            
            % Get f-measure
            for i=1:size(outputs,1)
                f_measure(inp,i,run,iter) = f_Measure(labels(testInd,:),outputs(i,:)');
            end
        end
    end
end
%%
median(median(f_measure,4),3);
%    0.8846    0.7562    0.8689
%    0.9505    0.9086    0.9421

fig(1) = figure(1);
boxplot(reshape(f_measure,10,18)');
grid on;
ax = gca;
ylabel('f-measure');
ax.TickLabelInterpreter = 'tex';
ax.XTickLabel = {'ANN_{p}','ANN_{s}','DT_{p}','DT_{s}','RF_{p}','RF_{s}','DP_{p}','DP_{s}','DM_{p}','DM_{s}'};
%title('Classification using parameter or similarity space');

fig(2) = figure(2);
for i=1:3
    subplot(3,1,i); hold off;
    boxplot(reshape(f_measure(:,:,:,i),10,6)'); hold on;
    grid on;ax = gca;ylabel('f-measure');ax.TickLabelInterpreter = 'tex';
    ax.XTickLabel = {'ANN_{p}','ANN_{s}','DT_{p}','DT_{s}','RF_{p}','RF_{s}','DP_{p}','DP_{s}','DM_{p}','DM_{s}'};
end

%title('Classification using parameter or similarity space');

if inclZero
    save_figures(fig, '.', 'classification_prototypes_incl0', 16, [6 4]);
else
    save_figures(fig, '.', 'classification_prototypes_no0', 16, [6 4]);
end



