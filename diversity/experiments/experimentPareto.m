load('neutrality400_DOF_16.mat');
%%

for r=1:length(radialBoundAdaptations)
    disp(r)
    for rep=1:replicates
        disp(rep)
        for pp=1:length(NSGA_polygons{r,rep})
            NSGA_polygons{r,rep}{pp}.Vertices = NSGA_polygons{r,rep}{pp}.Vertices-mean(NSGA_polygons{r,rep}{pp}.Vertices);
        end
        [NSGA_phenotypes{r,rep}] = getPhenotypeBoolean(NSGA_polygons{r,rep},resolution);
        for pp=1:length(QD_polygons{r,rep})
            QD_polygons{r,rep}{pp}.Vertices = QD_polygons{r,rep}{pp}.Vertices-mean(QD_polygons{r,rep}{pp}.Vertices);
        end
        [QD_phenotypes{r,rep}] = getPhenotypeBoolean(QD_polygons{r,rep},resolution);
        for pp=1:length(RLS_polygons{r,rep})
            RLS_polygons{r,rep}{pp}.Vertices = RLS_polygons{r,rep}{pp}.Vertices-mean(RLS_polygons{r,rep}{pp}.Vertices);
        end
        [RLS_phenotypes{r,rep}] = getPhenotypeBoolean(RLS_polygons{r,rep},resolution);
        
        
    end
end



%% Objective Space - Ground Truth Pareto Set

% Ground Truth Pareto Set
nShapes = 10
PARETO_genomes = [];
axDev = d.ranges(1,1):(range(d.ranges(1,:))./(nShapes-1)):d.ranges(1,2);
radDev = d.ranges(9,1):(range(d.ranges(9,:))./(nShapes-1)):d.ranges(9,2);
deletingIDs = find(abs(axDev)<0.1);
axDev(deletingIDs) = [];
radDev(deletingIDs) = [];
nShapes = length(axDev);
for i=1:nShapes
    for j=1:nShapes
        PARETO_genomes(end+1,:) = [axDev(i).*ones(1,DOF/2), radDev(j).*ones(1,DOF/2)];
    end
end

[PARETO_fitness,PARETO_polygons,~] = d.fitfun(PARETO_genomes);
PARETO_features = categorize(PARETO_polygons, d)';
[PARETO_phenotypes] = getPhenotypeBoolean(PARETO_polygons,resolution);

%% Distances
usePheno = true;
stdsNSGA = []; stdsQD = []; stdsRLS = []; stdsNSGA_t = []; stdsQD_t = []; stdsRLS_t = [];
if usePheno
    for r=1:length(radialBoundAdaptations)
        disp(['r: ' int2str(r)]);
        for rep=1:replicates
            disp(['rep: ' int2str(rep)]);
            for i=1:length(NSGA_polygons{r,rep})
                disp(['NSGA i: ' int2str(i)]);
                err = sum((xor(NSGA_phenotypes{r,rep}(i,:),PARETO_phenotypes))')./sum(PARETO_phenotypes');
                stdsNSGA_t(r,rep,i) = min(err);
            end
            for i=1:length(QD_polygons{r,rep})
                disp(['QD   i: ' int2str(i)]);
                err = sum((xor(QD_phenotypes{r,rep}(i,:),PARETO_phenotypes))')./sum(PARETO_phenotypes');
                stdsQD_t(r,rep,i) = min(err);
            end
            for i=1:length(RLS_polygons{r,rep})
                disp(['RLS  i: ' int2str(i)]);
                err = sum((xor(RLS_phenotypes{r,rep}(i,:),PARETO_phenotypes))')./sum(PARETO_phenotypes');
                stdsRLS_t(r,rep,i) = min(err);
            end
        end
    end
else
    for r=1:length(radialBoundAdaptations)
        disp(['r: ' int2str(r)]);
        for rep=1:replicates
            disp(['rep: ' int2str(rep)]);
            bla = pdist2(NSGA_genomes{r,rep},PARETO_genomes,'euclidean');
            stdsNSGA_t(r,rep,:) = min(bla')';
            bla = pdist2(QD_genomes{r,rep},PARETO_genomes,'euclidean');
            stdsQD_t(r,rep,:) = min(bla')';
            bla = pdist2(RLS_genomes{r,rep},PARETO_genomes,'euclidean');
            stdsRLS_t(r,rep,:) = min(bla')';
        end
    end
end


stdsNSGA = reshape(stdsNSGA_t,length(radialBoundAdaptations),[]);
stdsQD = reshape(stdsQD_t,length(radialBoundAdaptations),[]);
stdsRLS = reshape(stdsRLS_t,length(radialBoundAdaptations),[]);


%%
clear figs; close all;
set(0,'DefaultFigureWindowStyle','default')

allStds = 100*[stdsNSGA',stdsQD',stdsRLS'];
figs(1) = figure(1);hold off;
allStds = [reshape(allStds(:,1:5),1,[])',reshape(allStds(:,6:10),1,[])',reshape(allStds(:,11:15),1,[])'];
violins = violinplot(allStds, repelem([1:5],5*size(allStds,1)),'ViolinColor',[0 0 0],'ShowData',false)
grid on;
ax = gca;
%ax.YScale = 'log';
ax.XTick = 1:3;
ax.XTickLabel = {'NSGA-II','VE','RLS'}
axis([0 4, 0, 100]);
ylabel('% Pixel Error rel. to closest gt Pareto member');

%%
figs(end+1) = figure;
hold off;
nShapes = size(PARETO_genomes,1);
xPos = 0:ceil(sqrt(nShapes))-1; [X,Y] = ndgrid(xPos,xPos);
placement = [X(:) Y(:)]; placement = placement(1:nShapes,:);
showPhenotype(PARETO_genomes,d,gcf,placement);
axis equal;axis tight;
ax = gca;
ax.XTick = []; ax.YTick = [];
title('Ground Truth Pareto Set');

%%
cmap = redgreencmap(10+1,'Interpolation','linear');
r = 5
figs(end+1) = figure;
hold off;
ss = squeeze(stdsNSGA_t(r,rep,:));
[sortDist,sorted] = sort(ss,'ascend');
sortDist = floor(10*sortDist)+1;
nShapes = size(NSGA_genomes{r,rep},1);
xPos = 0:ceil(sqrt(nShapes))-1; [X,Y] = ndgrid(xPos,xPos);
placement = [X(:) ceil(sqrt(nShapes))-1-Y(:)]; placement = placement(1:nShapes,:);
figHandle = showPhenotype(NSGA_genomes{r,rep}(sorted,:),d,gcf,placement,cmap(sortDist,:));
axis equal;axis tight;
ax = gca;
ax.XTick = []; ax.YTick = [];
%title('NSGA-II, in order of closeness to Ground Truth');
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Pixel Error (%)';
cb.TickLabels = 0:10:100;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';


figs(end+1) = figure;
hold off;
ss = squeeze(stdsQD_t(r,rep,:));
[sortDist,sorted] = sort(ss,'ascend');
sortDist = floor(10*sortDist)+1;
nShapes = size(QD_genomes{r,rep},1);
xPos = 0:ceil(sqrt(nShapes))-1; [X,Y] = ndgrid(xPos,xPos);
placement = [X(:) ceil(sqrt(nShapes))-1-Y(:)]; placement = placement(1:nShapes,:);
showPhenotype(QD_genomes{r,rep}(sorted,:),d,gcf,placement,cmap(sortDist,:));
axis equal;axis tight;
ax = gca;
ax.XTick = []; ax.YTick = [];
%title('VE, in order of closeness to Ground Truth');
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Pixel Error (%)';
cb.TickLabels = 0:10:100;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';


figs(end+1) = figure;
hold off;
ss = squeeze(stdsRLS_t(r,rep,:));
[sortDist,sorted] = sort(ss,'ascend');
sortDist = floor(10*sortDist)+1;
nShapes = size(RLS_genomes{r,rep},1);
xPos = 0:ceil(sqrt(nShapes))-1; [X,Y] = ndgrid(xPos,xPos);
placement = [X(:) ceil(sqrt(nShapes))-1-Y(:)]; placement = placement(1:nShapes,:);
showPhenotype(RLS_genomes{r,rep}(sorted,:),d,gcf,placement,cmap(sortDist,:));
axis equal;axis tight;
ax = gca;
ax.XTick = []; ax.YTick = [];
%title('RLS, in order of closeness to Ground Truth');
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Pixel Error (%)';
cb.TickLabels = 0:10:100;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';


%%
save_figures(figs,'.','3algs_ParetoCloseness',12,[4 4])




%%
figure(1);
subplot(1,2,1);
imagesc(reshape(NSGA_phenotypes{r,rep}(2,:),64,64));
subplot(1,2,2);
imagesc(reshape(PARETO_phenotypes(2,:),64,64));

for i=1:10
100*sum(xor(NSGA_phenotypes{r,rep}(i,:),PARETO_phenotypes(i,:)))./sum(PARETO_phenotypes(i,:))
end

