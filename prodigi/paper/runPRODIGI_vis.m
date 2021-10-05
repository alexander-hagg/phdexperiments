%% Visualization
%% Latent Space
for it=1:length(output)
    disp(['t-SNE ' int2str(it) ]);
    if it == 1
        optima = reshape(output{it}.data{end}.sail.predMap(end).genes,625,10);
        optima(any(isnan(optima')),:) = [];
    else
        newOptima = reshape(output{it}.data{end}.sail.predMap(end).genes,625,10);
        newOptima(any(isnan(newOptima)'),:) = [];
        optima = [optima;newOptima];
    end
    wr
    [labels{it}, redX{it}, ~, ~, ~, clusters{it}] = dimReducedClustering( optima, 'tSNE', 2);
    oldLength{it} = size(optima,1);        
end

%%
set(0,'DefaultFigureWindowStyle','docked')
axisRanges = [-40 40 -40 40];
showIterations = 5;
dotsize = 12;
A = redX{1};
firstConceptIDs = ismember(output{1}.conceptLabels,output{1}.conceptSelection.id);

for i=1:showIterations
    
    B = redX{i};

    fig(i) = figure(i);hold off;
    scatter(B(size(redX{1},1)+1:end,1),B(size(redX{1},1)+1:end,2),dotsize/2,[0.5 0.5 0.5 ],'filled');
    hold on;
    scatter(B(1:size(redX{1},1),1),B(1:size(redX{1},1),2),dotsize,[1 0 0],'filled');
    selectedConceptCoords = B(1:size(redX{1},1),:);
    selectedConceptCoords = selectedConceptCoords(firstConceptIDs,:);
    scatter(selectedConceptCoords(:,1),selectedConceptCoords(:,2),dotsize,[0 1 0],'filled');
    
    legend('Optima', 'Non-Selected', 'Selected', 'Location','SouthEast');
    ax = gca;ax.XAxis.Visible = 'off';ax.YAxis.Visible = 'off';
    axis(axisRanges);
    title(['Optima/Categories']);
    
end
%save_figures(fig, './', ['prodigi_tsneOriginal_'], 16, [7 6]);

%% Correct by transformation with linear least squares on first set of optima
set(0,'DefaultFigureWindowStyle','default')

axisRanges = [-30 30 -40 40];
showIterations = 6;
dotsize = 48;
firstConceptIDs = ismember(output{1}.conceptLabels,output{1}.conceptSelection.id);

for i=1:showIterations
    
    B = redX{i};
    centroid_A = mean(A);
    centroid_B = mean(B(1:size(redX{1},1),:));
    N = size(A,1);
    H = (B(1:size(redX{1},1),:) - repmat(centroid_B, N, 1))' * (A - repmat(centroid_A, N, 1));
    [U,S,V] = svd(H);
    R = V*U';
    t = -R*centroid_A';% + centroid_B';
    B2 = (R*B') + repmat(t, 1, size(B,1));
    B2 = B2';
    
    fig(i) = figure(i);hold off;
    scatter(B2(size(redX{1},1)+1:end,1),B2(size(redX{1},1)+1:end,2),dotsize/4,[0.5 0.5 0.5],'.');hold on;
    scatter(B2(1:size(redX{1},1),1),B2(1:size(redX{1},1),2),dotsize,[0.5 0.5 0.5],'+');
    hold on;
    selectedConceptCoords = B2(1:size(redX{1},1),:);
    selectedConceptCoords = selectedConceptCoords(firstConceptIDs,:);
    scatter(selectedConceptCoords(:,1),selectedConceptCoords(:,2),dotsize,[0 0 0],'o','filled');
    
    %if i==1; legend('Non-Selected', 'Selected', 'Location','NorthWest');end
    if i==1; legend('New Optima', 'Non-Selected', 'Selected', 'Location','NorthWest');end
    ax = gca;ax.XAxis.Visible = 'off';ax.YAxis.Visible = 'off';
    axis(axisRanges);
    title(['Iteration ' int2str(i)]);
end

save_figures(fig, './', ['prodigi_tsneRotated_'], 18, [5 5]);


%% Show concept convex hulls
batchSizes = [0];
for l=1:length(labels)
    batchSizes(end+1) = length(labels{l});
end

for i=1:4%showIterations
    B = redX{i};
    centroid_A = mean(A);
    centroid_B = mean(B(1:size(redX{1},1),:));
    N = size(A,1);
    H = (B(1:size(redX{1},1),:) - repmat(centroid_B, N, 1))' * (A - repmat(centroid_A, N, 1));
    [U,S,V] = svd(H);
    R = V*U';
    t = -R*centroid_A' + centroid_B';
    B2 = (R*B');% + repmat(t, 1, size(B,1));
    B2 = B2';
    
    fig(i) = figure(i);hold off;
    concepts = unique(labels{i});
    for c=1:max(concepts)
        conceptMemberIDs = find(labels{i}==c);
        conceptMembers = B2(conceptMemberIDs,:);
        if size(conceptMembers,1)>2
            K = boundary(conceptMembers(:,1),conceptMembers(:,2));
            fill(conceptMembers(K,1),conceptMembers(K,2),'r','LineWidth',1); hold on;
        end
    end
    
    ax = gca;ax.XAxis.Visible = 'off';ax.YAxis.Visible = 'off';
    axis(axisRanges);
end

%save_figures(fig, './', ['prodigi_tsnePrototypes_'], 16, [7 6]);





%% View Maps (Predicted and True Fitness, Fitness improvement, Model Eror)
set(0,'DefaultFigureWindowStyle','docked')

for it=1:length(output)
    figure(1);
    subplot(p.numIterations,4,(it-1)*4+1);
    viewMap(output{it}.data{end}.sail.predMap(end).fitness,d);title('Predicted Fitness');caxis([-5 -3]);
    subplot(p.numIterations,4,(it-1)*4+2);
    viewMap(output{it}.data{end}.sail.predMap(end).fitness_true,d);title('True Fitness');caxis([-5 -3]);
    if it > 1
        subplot(p.numIterations,4,(it-1)*4+3);
        viewMap(output{selectedRUN}.data{end}.sail.predMap(end).fitness_true-output{it}.data{end}.sail.predMap(end).fitness_true,d);title('True Fitness Improvement');caxis([-0.5 0.5]);
    end
    subplot(p.numIterations,4,(it-1)*4+4);
    err = abs(output{it}.data{end}.sail.predMap(end).fitness-output{it}.data{end}.sail.predMap(end).fitness_true);
    viewMap(err,d);title('Abs. Error');caxis([0 1]);
end

%% Boxplots (True Fitness, Improvement, Abs. Model Error)
set(0,'DefaultFigureWindowStyle','default')
for it=1:length(output)
    fits(it,:) = reshape(output{it}.data{end}.sail.predMap(end).fitness,625,1);
    fittrues(it,:) = reshape(output{it}.data{end}.sail.predMap(end).fitness_true,625,1);
    if it > 1, improv(it,:) = fittrues(it-1,:)-fittrues(it,:);end
    errs(it,:) = abs(fits(it,:)-fittrues(it,:));
end

figure(1);
boxplot(fittrues'); axis([0.5 [length(output) + 0.5] -5.2 -3.5]);
title('True Fitness');grid on;

figure(2);
boxplot(improv'); axis([0.5 [length(output) + 0.5] -0.5 1]);grid on;
title('True Fitness Improvement');

figure(3);
boxplot(errs'); axis([0.5 [length(output) + 0.5] -0.01 0.7]);grid on;
title('Abs. Error');

