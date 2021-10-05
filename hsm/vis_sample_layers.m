function [ output_args ] = vis_sample_layers( model,layer_prediction,prediction_ids,inputs,outputs )
%% VIS_SAMPLE_LAYERS Visualizes a hierarchical model's sample subsets per layer


%% model.hierarchy.Node{m}{2}
symbols = {'x','o','+','*','s','d','<','>','x','o','+','*','s','d','<','>','x','o','+','*','s','d','<','>','x','o','+','*','s','d','<','>'};
layers = get_layers(model);

%%
figure(1);
for l=1:length(layers)
    subplot(1,length(layers),l);hold off;
    allsamples = ones(size(outputs));
    for m=1:length(layers{l})
        subset = prediction_ids{l} == layers{l}(m);
        allsamples = allsamples - subset;
        plot(outputs(subset),layer_prediction{l}(:,subset),symbols{m});
        hold on;
    end
    if sum(allsamples) > 0;plot(outputs(:,logical(allsamples)),layer_prediction{l}(:,logical(allsamples)),'x');hold on;end;
    mdl = fitlm(layer_prediction{l},outputs);
    
    plot([-1.5 1.5],[-1.5 1.5]);hold on;
    grid on; axis equal;
    axis([-1.2 1.2 -1.2 1.2]);
    axis([-1.5 1.5 -1.5 1.5]);
    
    xlabel('Target');ylabel('Output');
    title(['Layer ' int2str(l) ', R^2: ' num2str(mdl.Rsquared.Ordinary)]);
    hold off;
end


%%

figure(2);
for l=1:length(layers)
    subplot(1,length(layers),l);hold off;
    allsamples = ones(size(outputs));
    plot(inputs(1,:),outputs,'-x');hold on;
    for m=1:length(layers{l})
        subset = prediction_ids{l} == layers{l}(m);
        allsamples = allsamples - subset;
        plot(inputs(1,subset),layer_prediction{l}(:,subset),symbols{m});
    end
    if sum(allsamples) > 0;plot(inputs(:,logical(allsamples)),layer_prediction{l}(:,logical(allsamples)),'x');hold on;end;
    
    grid on; axis equal;
    axis([-1.2 1.2 -1.2 1.2]);
    axis([-1.5 1.5 -1.5 1.5]);
    
    xlabel('Target');ylabel('Output');
    title(['Predictions Layer ' int2str(l)]);
end
hold off;



end

