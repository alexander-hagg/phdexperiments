function fig = vPlans(genomes,d,varargin)
%VPLANS Summary of this function goes here
%   fig = vPlans(genomes,d,varargin)
%
d.map = importdata(['domains/mazeDomain/simulator/worlds/' d.mazeCfgBaseFile '.pbm']);

if nargin > 2
    if ~isempty(varargin{1})
        fig = figure(varargin{1});
    else
        fig = figure;
    end
else
    fig = figure;
end

if nargin > 3
    clrs = varargin{2};
else
    clrs = repmat([0 0 0],size(genomes,1),1);
end

if isempty(fig.Children)
    hold off;
    disp('draw background');
    bg = imagesc(d.map); hold on; colormap(colorcube(2));
end



%%
validInds = true(1, size(genomes,1));
% Calculate validity
if strcmp(d.representation,'planner');
    [validInds,phenotypes] = feval(d.validate,genomes,d);
    clrsVal = clrs(validInds,:);
    valPhen = phenotypes(validInds,:,:);
    for i=1:size(valPhen,1)
        plot(squeeze(valPhen(i,1,:))',squeeze(valPhen(i,2,:))','LineWidth',1,'Color',clrsVal(i,:));
        hold on;
    end
    
    plot(squeeze(phenotypes(~validInds,1,:))',squeeze(phenotypes(~validInds,2,:))','LineWidth',8,'Color',[0.2 0.2 0.2]);
    
    xx = squeeze(phenotypes(:,1,[1,end]));
    yy = squeeze(phenotypes(:,2,[1,end]));
    
    scatter(xx(:),yy(:),32,[repmat([0 0 0],length(validInds),1);clrs],'filled');
    
else
    [phenotypes,~] = eval_maze(genomes,d);
    
    clrsVal = clrs(validInds,:);
    valPhen = phenotypes(validInds,:,:);
    for i=1:size(valPhen,1)
        plot(squeeze(valPhen(i,:,1))',squeeze(valPhen(i,:,2))','LineWidth',1,'Color',clrsVal(i,:));
        hold on;
    end
    
    plot(squeeze(phenotypes(~validInds,:,1))',squeeze(phenotypes(~validInds,:,2))','LineWidth',8,'Color',[0.2 0.2 0.2]);
    
    xx = squeeze(phenotypes(:,[1,end],1));
    yy = squeeze(phenotypes(:,[1,end],2));
    
    scatter(xx(:),yy(:),32,[repmat([0 0 0],length(validInds),1);clrs],'filled');
    
end

axis([0 400 -0 400]);

%% Draw selection circles
viscircles(d.goalRings{1},[30 30 30]); hold on;
viscircles(d.center,d.ringDiameter(1)/2); hold on;


end

