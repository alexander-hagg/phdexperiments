clear

% Configure experiment
visOnline = false;
goals = [300,10;490,500;110,500];
numWeights = 23;
timesteps = 150;
maze = 'hardRound'; %'exampleRound'
fitfun = @(data) min(sqrt(sum(((goals-data(end,1:2)).^2)')));


% Configure GA
numGen = 256;
popsize = 32;
mutDist = 0.5;
mutProb = 1/numWeights;
weights = 2*rand(popsize,numWeights)-1;

% Configure dimensionality reduction
mapMethod = 'SNE'; %'tSNE' 'AutoEncoder'


for gen=1:numGen
    % Evaluate
    [fitness,trajectories] = eval_maze(weights,maze,timesteps,fitfun,false);
    
    % Elitism
    allFit(gen,:) = fitness';
    [eliteFIT(gen),eliteID] = min(fitness);
    disp(['Current fit: ' num2str(eliteFIT(gen))]);
    
    % Tournament selection
    tourIDs = randi(popsize,2,2*popsize);
    tours = fitness(tourIDs);
    winnerIDs =  [~(tours(2,:)<=tours(1,:));(tours(2,:)<=tours(1,:))];
    winners = tourIDs(winnerIDs);
    
    % Recombination
    % Crossover
    w1 = weights(winners(1:end/2),:);
    w2 = weights(winners(end/2+1:end),:);
    unif = logical(randi(2,size(w1,1),size(weights,2))-1);
    newWeights = zeros(size(weights));
    newWeights(unif) = w1(unif);
    newWeights(~unif) = w2(~unif);
    % Mutation
    doMut = rand(size(newWeights)) < mutProb;
    newWeights = newWeights + mutDist*randn(size(newWeights)).*doMut;
    
    % Rest elitism
    newWeights(1,:) = weights(eliteID,:);
    weights = newWeights;
    
    %% Visualization
    if visOnline
        figure(1); hold off;
        plot(eliteFIT,'LineWidth',2);
        hold on;
        plot(median(allFit'),'r','LineWidth',2);
        legend('Elite','Median');
        grid on;
        if gen > 1; axis([1 length(eliteFIT) 0 300]);end
        xlabel('Gens');ylabel('Distance to *a* goal');
        
        [reducedParSpace, mapping] = compute_mapping([weights;zeros(1,size(weights,2));ones(1,size(weights,2));-ones(1,size(weights,2))], mapMethod, 2);
        figure(2); hold off;
        colors = parula(4);
        sc1 = scatter(reducedParSpace(1:end-3,1),reducedParSpace(1:end-3,2),64,colors([ones(1,size(weights,1))]),'filled');
        hold on;
        sc2 = scatter(reducedParSpace(end-2:end,1),reducedParSpace(end-2:end,2),64,colors([2 2 2]),'filled');
        title(['Weights - ' mapMethod]);
        legend([sc1,sc2],'Weights', '[-1,...,-1] [0,...,0] [1,...,1]');
        drawnow;
    end
end

%%

[fitness,trajectories] = eval_maze(weights,maze,timesteps,fitfun,true);
%%
fig = visTrajectories(maze, trajectories, 99);



%%
colors = parula(4);
allweights = [weights1;weights2;weights3;weights4];
colorIDs = [ones(1,size(weights1,1)), 2*ones(1,size(weights2,1)), 3*ones(1,size(weights3,1)), 4*ones(1,size(weights4,1))];

mapMethod = 'SNE'; %'tSNE' 'AutoEncoder'
figure(3);
subplot(3,1,1);
hold off;
[reducedParSpace, mapping] = compute_mapping(allweights, mapMethod, 2);
scatter(reducedParSpace(:,1),reducedParSpace(:,2),64,colors(colorIDs,:),'filled');
title(['Weights of 4 runs: ' mapMethod]);
hold on;
%
mapMethod = 'tSNE'; %'tSNE' 'AutoEncoder'
figure(3);
subplot(3,1,2);
hold off;
[reducedParSpace, mapping] = compute_mapping(allweights, mapMethod, 2, size(allweights, 2), 50);
scatter(reducedParSpace(:,1),reducedParSpace(:,2),64,colors(colorIDs,:),'filled');
title(mapMethod);
hold on;
%
mapMethod = 'AutoEncoder'; %'tSNE' 'AutoEncoder'
figure(3);
subplot(3,1,3);
hold off;
[reducedParSpace, mapping] = compute_mapping(allweights, mapMethod, 2);
scatter(reducedParSpace(:,1),reducedParSpace(:,2),64,colors(colorIDs,:),'filled');
title(mapMethod);
hold on;
%axis([-1000 1000 -1000 1000]);

