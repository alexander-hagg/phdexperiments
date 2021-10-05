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

clear deviation classMembership;
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
        
        % Train GP mapping
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
        outputs = [m_1,m_2];
        conf_outputs = [s2_1,s2_2];
        
        % Get sets of candidate classes (0.5, 1, 1.5, 2, 2.5, 3 sigma threshold)
        sigmas = [0.0, 0.5,1,1.5,2,2.5,3];
        for sigID = 1:length(sigmas)
            % Get sets of class members
            lbls = output{run}.data{iter}.conceptLabels;
            clear deltaClasses;
            for i=1:max(lbls)
                memberIDs = (lbls==i);
                members = latent(memberIDs,:);
                % Get distance to class for all outputs
                deltaClasses(i,:) = min(pdist2(outputs,members)');
            end
            delta = deltaClasses - min(deltaClasses); % Closest class at zero
            variance = sqrt(max(conf_outputs'));
            classMembership{run,iter}(sigID,:,:) = bsxfun(@le,delta,sigmas(sigID)*variance);
        end
    end
end

%% Count number of classes
countClasses = [];
for run=2:length(output)
    for iter=1:length(output{run}.data)
        countClasses = [countClasses,squeeze(sum(classMembership{run,iter},2))];
    end
end

fig(1) = figure(1);hold off;
boxplot(countClasses','Symbol','k');
grid on;
ax = gca;
ax.XTickLabels = {'0\sigma', '0.5\sigma','1\sigma','1.5\sigma','2\sigma','2.5\sigma','3\sigma'};
ylabel('\# Class Assignments');
xlabel('Constraint Threshold');
axis([0.5 7.5 -0.2 8]);
set(gca,'TickLabelInterpreter','tex');
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
hold on;
sc = scatter([1 2 3 4 5 6 7],mean(countClasses'),'filled');
legend(sc,'Mean');

save_figures(fig, '.', 'softcontraint_boxplots', 16, [6 4]);
