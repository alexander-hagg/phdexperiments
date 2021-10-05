thetaValues = [30,90,150,210,270,330];
dimReductionMethods = {'none','PCA','Isomap','KernelPCA','SNE', 'tSNE','AutoEncoder'};
%dimReductionMethods = {'none','PCA','tSNE'};

for thetaID=1:length(thetaValues)
    disp(['Theta: ' int2str(thetaValues(thetaID))]);
    
    load(['PROPHESAI/ECJ_THETA_baselines/HID5_64kgens/' int2str(thetaValues(thetaID)) '.mat']);
    genes = runData{1}.acqMap.genes;
    genes = reshape(genes,size(genes,1)*size(genes,2),size(genes,3));
    genes(any(isnan(genes')),:) = [];
    
    for dr=7:length(dimReductionMethods)
        disp(['DR method: ' dimReductionMethods{dr}]);
        [classification{dr,thetaID},stats{dr,thetaID}] = extractClasses(genes,dimReductionMethods{dr},10);
    end
    
end

save('data_analysis_dimensionalityreduction_neural.mat');
%% Analysis
for dr=1:length(dimReductionMethods)
    valGPLUS{dr} = []; valGPLUS_ORG{dr} = [];
    for thetaID=1:length(thetaValues)
        %classification{dr,thetaID}
        valGPLUS{dr} = [valGPLUS{dr}  stats{dr,thetaID}.valGPLUS];
        valGPLUS_ORG{dr} = [valGPLUS_ORG{dr}  stats{dr,thetaID}.valGPLUS_ORG];
        
    end
end


%% Visualization
dimReductionMethods = {'none','PCA','Iso','kPCA','SNE', 'tSNE','AutoEnc'};

fig(1) = figure(1);
boxplot(log([valGPLUS{1}' valGPLUS{2}' valGPLUS{3}' valGPLUS{4}' valGPLUS{5}' valGPLUS{6}' valGPLUS{7}']),'PlotStyle','compact')
axis([0.5 7.5 -16 -7]);
ax = gca;
ax.XTick = 1:7
ax.XTickLabel = dimReductionMethods;
grid on; grid minor;
ylabel('log(G+)');
ax.XTickLabelRotation = 90

fig(2) = figure(2);
boxplot(log([valGPLUS_ORG{1}' valGPLUS_ORG{2}' valGPLUS_ORG{3}' valGPLUS_ORG{4}' valGPLUS_ORG{5}' valGPLUS_ORG{6}' valGPLUS_ORG{7}']),'PlotStyle','compact')
axis([0.5 7.5 -16 -7]);
ax = gca;
ax.XTick = 1:7
ax.XTickLabel = dimReductionMethods;
grid on; grid minor;
ylabel('log(G+)');
ax.XTickLabelRotation = 90

save_figures(fig, '.', 'drgplus', 12, [3 4])

%%

theta=4;dr = 6;
maxLabel = max(classification{dr,theta}.labels);
colors = hsv(maxLabel);
scatter(classification{dr,theta}.simX(:,1),classification{dr,theta}.simX(:,2),32,colors(classification{dr,theta}.labels,:),'filled')
title([int2str(maxLabel) ' classes']);
grid on; grid minor;

