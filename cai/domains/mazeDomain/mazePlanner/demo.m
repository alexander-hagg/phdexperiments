clear;

% Genome: 2 X N, N path components
% angle, length
numSamples = 30;
pathLength = 10;
genome = rand(numSamples,pathLength*2);
start = [200,200]';
map = importdata('mediumRound.pbm');
d.ranges                    = [-pi -20; pi 20];

d.angle = 30;
d.angle = (d.angle)/180 * pi;


tic;
% Calculate phenotype
for i=1:numSamples
    phenotype(i,:,1) = start;
    heading = d.angle;
    for p=1:pathLength
        a = genome(i,p) * range(d.ranges(:,1)) - d.ranges(2,1);
        l = genome(i,p+pathLength) * range(d.ranges(:,2)) - d.ranges(2,2);
        heading = heading + a;
        
        dx = l / sin(heading);
        dy = l / cos(heading);
        
        phenotype(i,:,p+1) = phenotype(i,:,p) + [dx;dy]';
    end
    
    % Calculate validity
    [cx{i},cy{i},c{i}] = improfile(map,phenotype(i,1,:),phenotype(i,2,:));
    cx{i} = cx{i}(~isnan(c{i}));
    cy{i} = cy{i}(~isnan(c{i}));
    c{i} = c{i}(~isnan(c{i}));
end

toc

%%

figure(1);hold off;
bg = imagesc(map); hold on; colormap(colorcube(2));
for i=1:length(cx)
    disp(i)
    for p=1:length(cx{i})-1
        pl = plot(cx{i}(p:p+1),cy{i}(p:p+1),'LineWidth',4,'Color',[0 1/i i/5].*c{i}(p)+[1 0 0].*(1-c{i}(p)));
        if p==1; plots(i) = pl;end
    end
end
axis([0 400 -0 400]);
legend(plots,'1','2');
drawnow;

%%
bg = imagesc(d.map); hold on; colormap(colorcube(2));
for i=1:size(phenotype,1)
    plot(squeeze(phenotype(i,:,1)),squeeze(phenotype(i,:,2)));
    axis([0 400 -0 400]);
end
    drawnow;





%%



