numHiddenS = [2 5 10];
for nH=1:length(numHiddenS)
    numHidden=numHiddenS(nH);
    numThreshs = [0.9 0.1];
    for threshID=1:length(numThreshs)
        thresh = numThreshs(threshID);
        load(['out_' int2str(numHidden) '_thresh_' num2str(thresh) '_iteration_5.mat']);
        for iter=1:5
            alignedVals = [runData{iter}.classification.values{1,:}]; 
            popsize(nH,threshID,iter) = size(alignedVals,2);
            stats(nH,threshID,iter) = sum(alignedVals<1)/size(alignedVals,2);
        end
    end
end

%%
for nH=1:length(numHiddenS)
    figure(nH);plot(100*squeeze(stats(nH,:,:))','LineWidth',2);
    grid on;grid minor;
    axis([1 5 50 90]);
    ylabel('% forward moving');xlabel('Class Selection at Generation');
    ax = gca;xticks(1:5);xticklabels({'4000','5000','6000','7000','8000'});
    title(['Num Hidden Neurons: ' int2str(numHiddenS(nH))]);
    legend('Select classes >90% fwd moving','Select classes <10% fwd moving','Location','NorthWest');
    
    figure(nH+length(numHiddenS));plot(squeeze(popsize(nH,:,:))','LineWidth',2);
    grid on;grid minor;
    axis([1 5 300 500]);
    ylabel('# Solutions Found');xlabel('Class Selection at Generation');
    ax = gca;xticks(1:5);xticklabels({'4000','5000','6000','7000','8000'});
    title(['Num Hidden Neurons: ' int2str(numHiddenS(nH))]);
    legend('Select classes >90% fwd moving','Select classes <10% fwd moving','Location','NorthWest');
        
end