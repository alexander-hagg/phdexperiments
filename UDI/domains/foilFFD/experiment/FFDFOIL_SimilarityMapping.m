% Similarity Mapping, (2)
% Hard Constraints (1)
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

%inclZero = false;


clear deviation;
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 0.15;

for run=2:length(output)
    disp(['Run ' int2str(run)]);
    for iter=1:length(output{run}.data)
        %%
        disp(['Iter ' int2str(iter)]);
        % Get inputs
        optima = output{run}.data{iter}.optima;
        latent = output{run}.data{iter}.latent;
        Q = size(optima,1);
        % Train test split
        [trainInd,valInd,testInd] = dividerand(Q,trainRatio,valRatio,testRatio);
        
        clear outputs conf_outputs;
        
        tic1 = tic;
        % ANN
        % Create a Pattern Recognition Network
        hiddenLayerSize = 25;
        net = feedforwardnet([hiddenLayerSize,hiddenLayerSize]); net.trainParam.showWindow = 0;
        % Set up Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = trainRatio + valRatio;net.divideParam.valRatio = valRatio;net.divideParam.testRatio = 0;
        % Train the Network
        [net,tr] = train(net,optima([trainInd,valInd],:)',latent([trainInd,valInd],:)');
        % Test the Network
        outputs(1,:,:) = net(optima(testInd,:)')';
        timings{run,iter,1} = toc(tic1)
        
        tic2 = tic;
        % GP
        p.covfunc   = @covSEard;             
        p.hyp.cov   = [zeros(size(optima,2),1);0]; % (unit vector in log space)
        p.meanfunc  = {@meanConst};  
        p.hyp.mean  = 0;
        p.likfunc   = @likGauss;     
        p.hyp.lik   = log(0.1);
        p.functionEvals = 100;      % function evals to optimize hyperparams
        GP_model1.hyp = minimize_gpml(p.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, ...
            p.covfunc, p.likfunc, optima([trainInd,valInd],:), latent([trainInd,valInd],1));
        GP_model2.hyp = minimize_gpml(p.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, ...
            p.covfunc, p.likfunc, optima([trainInd,valInd],:), latent([trainInd,valInd],2));        
        trainInput  = optima([trainInd,valInd],:);
        trainOutput = latent([trainInd,valInd],:);
        [m_1, s2_1] = gp(GP_model1.hyp, @infExact, p.meanfunc, p.covfunc, p.likfunc, ...
            trainInput, trainOutput(:,1), optima([testInd],:));
        [m_2, s2_2] = gp(GP_model2.hyp, @infExact, p.meanfunc, p.covfunc, p.likfunc, ...
            trainInput, trainOutput(:,2), optima([testInd],:));
        outputs(2,:,:) = [m_1,m_2];
        conf_outputs(1,:,:) = [s2_1,s2_2];        
        timings{run,iter,2} = toc(tic2)

        tic3 = tic;
        % GP
        p.covfunc   = {@covMaternard, 3};
        p.hyp.cov   = [zeros(size(optima,2),1);0]; % (unit vector in log space)
        p.meanfunc  = {@meanConst};  
        p.hyp.mean  = 0;
        p.likfunc   = @likGauss;     
        p.hyp.lik   = log(0.1);
        p.functionEvals = 100;      % function evals to optimize hyperparams
        GP_model1.hyp = minimize_gpml(p.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, ...
            p.covfunc, p.likfunc, optima([trainInd,valInd],:), latent([trainInd,valInd],1));
        GP_model2.hyp = minimize_gpml(p.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, ...
            p.covfunc, p.likfunc, optima([trainInd,valInd],:), latent([trainInd,valInd],2));        
        trainInput  = optima([trainInd,valInd],:);
        trainOutput = latent([trainInd,valInd],:);
        [m_1, s2_1] = gp(GP_model1.hyp, @infExact, p.meanfunc, p.covfunc, p.likfunc, ...
            trainInput, trainOutput(:,1), optima([testInd],:));
        [m_2, s2_2] = gp(GP_model2.hyp, @infExact, p.meanfunc, p.covfunc, p.likfunc, ...
            trainInput, trainOutput(:,2), optima([testInd],:));
        outputs(3,:,:) = [m_1,m_2];
        conf_outputs(2,:,:) = [s2_1,s2_2];
        timings{run,iter,3} = toc(tic3)   
        
        
        for i=1:size(outputs,1)
            %immse(latent(testInd,:),squeeze(outputs(i,:,:))).^(-0.5)
            out = abs(sum(bsxfun(@minus,latent(testInd,:),squeeze(outputs(i,:,:))),2))./range(latent(:));
            allOutputs{run,iter,i} = squeeze(outputs(i,:,:));
            if i>1
                allConfidences{run,iter,i-1} = squeeze(conf_outputs(i-1,:,:));        
            end
            gt{run,iter,i} = latent([testInd],:);
            if run==2 && iter == 1
                deviation{i} = out;
            else
                deviation{i} = [deviation{i};out];
            end
        end
        
    end
end
deviation = cell2mat(deviation);

%%

fig(1) = figure(1);
boxplot(100*deviation,'Symbol','k');
grid on;
ax = gca;
ax.XTickLabels = {'ANN', 'GP-SE-ard', 'GP-Matern-ard'};
ylabel('Rel. Abs. Error [\%]');
axis([0.5 3.5 -0.2 21]);
set(gca,'TickLabelInterpreter','latex');
