fig(1) = figure(99);
    
for i=1:length(output)
    times(i) = cell2mat(output{i}.trainTime(1,:));
end
times = times./p.nInitialSamples;

plot(numDims,times','o-');        
xlabel('Dimensions');
ylabel('Seconds');
title('Training Times');
grid on;


%%
fig(2) = figure(100);
    
for i=1:length(output)
    times(i) = cell2mat(output{i}.predictTime(1,:));
end
times = times./p.nInitialSamples;
plot(numDims,times','o-');        
xlabel('Dimensions');
ylabel('Seconds');
title('Prediction Times');
grid on;


%% Save figures

save_figures(fig, '.', 'timings', 18)

