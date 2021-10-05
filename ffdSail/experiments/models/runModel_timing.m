clear;clc;
domainname = 'PARSEC';
cd('../..');systemInit;

%%
if strcmp(domainname,'PARSEC');d = parsec_Domain; end
if strcmp(domainname,'FOILFFD');d = ffd_Domain; end
[~,~,~,modelCfgs] = models_LoadCfg(d);

p = sail;
s_numSamples       = [10 50 100 500 1000 2000];
stats.time.trainingSamples = s_numSamples;
        
for sSize=1:length(s_numSamples)
    numSamples = s_numSamples(sSize);
    %[observation, value] = initialSampling(d,numSamples);
    
    xp(sSize).gens = 500;
    xp(sSize).children = 100;
    xp(sSize).finalSamples = numSamples;
    xp(sSize).samplesAddedInIteration = 10;
    testObservation = rand(xp(sSize).gens*xp(sSize).children*(xp(sSize).finalSamples/xp(sSize).samplesAddedInIteration),10);
    stats.time.predictionSamples(sSize) = size(testObservation,1);

    disp(['Testing model ' modelCfgs{1}{1} ' with ' int2str(numSamples) ' training samples and ' int2str(size(testObservation,1)) ' test samples']);
    
    % Set acquisition and prediction model parameters
    for target=1:2    
        d.params{target}        = feval(['params' modelCfgs{1}{2}], modelCfgs{1}{3:end});
    end

    for target=1:2
        runTime = tic;
        model{target} = feval(['train' d.params{target}.name], observation, value(:,target), d.params{target});    
        stats.time.training(sSize, target) = toc(runTime);
        
        runTime = tic;        
        feval(['predict' d.params{target}.name], model{target}, testObservation);
        stats.time.prediction(sSize, target) = toc(runTime);
        
    end
end
save([resultpath '/results_MODELTIMING_' domainname '_' modelCfgs{1}{1} '.mat'], 'stats', 'p', 'd', 'xp', 'model', '-v7.3');

%%
fig(1) = figure(1); hold off;
loglog(stats.time.trainingSamples, mean(stats.time.training,2), 'r', 'LineWidth',2);
hold on;
ax1 = gca;
ax1.XTick = stats.time.trainingSamples;  
loglog(stats.time.trainingSamples, mean(stats.time.prediction,2), 'LineWidth',2, 'Color', 'k');

grid on;
ylabel('time [s]');
xlabel('real evaluations');

ax = gca;
ax.XTick = stats.time.trainingSamples;
axis([10 2000 1e-2 1e3]);
title('Total Training and Prediction Time');
legend('Training', 'Prediction', 'Location', 'SouthEast');

save_figures(fig, './', ['GP_times_'], 12, [7 5]);
