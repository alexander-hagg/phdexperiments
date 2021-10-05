clear
set(0,'DefaultFigureWindowStyle','normal')
rep = {'Parsec','FFD','CPPN'};
for iRep=1:3;load([rep{iRep} '_result']);end

%--- For plotting only ---|
p = sail; d = ffd_Domain; 
load(['~/Code/data/ffdSail/sailCPPN_true.mat']);
edges = output{1}.predMap(end).edges;
evalPt= Parsec.evalPt;
%-------------------------|

byItr = @(x) reshape(x,[625 20]);
colors = parula(4);
%% Accuracy of Models
fh(1) = figure(1); clf; fh1.Name = 'Model Accuracy'; hold on;
fh(1).Position = [2751,151,850,902];
% at each PE
acc{1} = nan(625,80); acc{1}(:,1:4:80) = byItr(Parsec.absErrorcD.^2);
acc{2} = nan(625,80); acc{2}(:,2:4:80) = byItr(FFD.absErrorcD.^2);
acc{3} = nan(625,80); acc{3}(:,3:4:80) = byItr(CPPN.absErrorcD.^2);

for i=1:3
    h = boxplot(acc{i},'PlotStyle','compact','Color',colors(i,:),'Labels',repmat({' '},[1,80]),'whisker',2);
end
xticks([2:4:80]);xticklabels(evalPt); set(gca,'YScale','log','YLim',[1e-13 1]); grid on;
xlabel('Precise Evaluations'); ylabel('MSE');

for i = 1:3; plot(NaN,1,'color', colors(i,:), 'LineWidth', 4); hold on;end;
legend(rep,'Location','SouthWest','FontSize',16);


title('Model Performance Per PE','FontSize',20); hold off;

%% at 1000 PE (Maps)
fh(2) = figure(2); clf; fh2.Name = 'Model Accuracy at 1000 PE'; hold on;
fh(2).Position = [1923,151,1221,943];
dragRange = [-5 -3.5]; absErrorRange = [0 0.5]; perErrorRange = [0 15];
subplot(3,3,1); 
    viewMap(Parsec.predcD(:,:,end),d,edges,'flip'); 
    caxis(dragRange);
    title('Parsec Drag Prediction')
subplot(3,3,2);
    viewMap(abs(Parsec.absErrorcD(:,:,end)),d,edges,'flip'); 
    caxis(absErrorRange);
    title('Parsec Absolute Drag Error')
subplot(3,3,3);
    viewMap(Parsec.percErrorcD(:,:,end),d,edges,'flip'); 
    caxis(perErrorRange);
    title('Parsec Percentage Drag Error')
    
subplot(3,3,4); 
    viewMap(FFD.predcD(:,:,end),d,edges); caxis(dragRange)
    title('FFD Drag Prediction')
subplot(3,3,5);
    viewMap(abs(FFD.absErrorcD(:,:,end)),d,edges); caxis(absErrorRange)
    title('FFD Absolute Drag Error')
subplot(3,3,6);
    viewMap(FFD.percErrorcD(:,:,end),d,edges); caxis(perErrorRange)
    title('FFD Percentage Drag Error')
    
subplot(3,3,7); 
    viewMap(CPPN.predcD(:,:,end),d,edges); caxis(dragRange)
    title('CPPN Drag Prediction')
subplot(3,3,8);
    viewMap(abs(CPPN.absErrorcD(:,:,end)),d,edges); caxis(absErrorRange)
    title('CPPN Absolute Drag Error')
subplot(3,3,9);
    viewMap(CPPN.percErrorcD(:,:,end),d,edges); caxis(perErrorRange)
    title('CPPN Percentage Drag Error')


%% Performance
%% at each PE
stdDevScale = 4;

fh(3) = figure(3); clf; fh3.Name = 'Map Performance per PE'; hold on;
fh(3).Position = [2751,151,850,902];

% Legend
for i = 1:3; plot(NaN,1,'color', colors(i,:), 'LineWidth', 4); hold on;end;
legend(rep,'Location','SouthWest','FontSize',16);

% Data
[l,p] = boundedline(evalPt,Parsec.trueFitOptMedian,Parsec.trueFitOptStd./stdDevScale);
l.Color = colors(1,:); l.LineWidth = 4; p.FaceColor = colors(1,:); p.FaceAlpha = 0.25;

[l,p] = boundedline(evalPt,FFD.trueFitOptMedian,FFD.trueFitOptStd./stdDevScale);
l.Color = colors(2,:); l.LineWidth = 4; p.FaceColor = colors(2,:); p.FaceAlpha = 0.25;

[l,p] = boundedline(evalPt,CPPN.trueFitOptMedian,CPPN.trueFitOptStd./stdDevScale);
l.Color = colors(3,:); l.LineWidth = 4; p.FaceColor = colors(3,:); p.FaceAlpha = 0.25;

grid on; set(gca, 'XLim',[100 1000])
title('Median Map Performance Per PE','FontSize',20); hold off;
xlabel('Precise Evaluations')
ylabel('Percentage of Optima')


%% at 1000 PE (Maps)
fh(4) = figure(4); clf; fh4.Name = 'Map Performance at 1000 PE'; hold on;
fh(4).Position = [1923,151,767,943];

%dragRange = [-5 -2.5];
dragRange = [-5 0];
subplot(3,2,1); 
    viewMap(Parsec.predFitness(:,:,end),d,edges,'flip'); 
    caxis(dragRange); 
    title('Parsec Fitness Prediction')
subplot(3,2,2);
    viewMap(Parsec.trueFitness(:,:,end),d,edges,'flip'); 
    caxis(dragRange);
    title('Parsec True Fitness')
    
subplot(3,2,3); 
    viewMap(FFD.predFitness(:,:,end),d,edges); 
    caxis(dragRange)
    title('FFD Fitness Prediction')
subplot(3,2,4);
    viewMap(FFD.trueFitness(:,:,end),d,edges); 
    caxis(dragRange)
    title('FFD True Fitness')
    
subplot(3,2,5); 
    viewMap(CPPN.predFitness(:,:,end),d,edges); 
    caxis(dragRange)
    title('CPPN Fitness Prediction')
subplot(3,2,6);
    viewMap(CPPN.trueFitness(:,:,end),d,edges); 
    caxis(dragRange)
    title('CPPN True Fitness')

colormap(parula(15));


%% Selection of Designs

% Maps

% Designs

%% Save all figures
figure(1); save2pdf('modelAccuracy');
figure(2); save2pdf('modelMaps');
figure(3); save2pdf('performance');
figure(4); save2pdf('perfMaps');
    