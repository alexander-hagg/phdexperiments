d = mirror_Domain; % for point IDs
numSamples = 100000;

samples = rand(numSamples,51);
tic;
for i=1:numSamples
    FV{i} = mirror_ffd_Express(samples(i,:), d.FfdP);
end
time(1) = toc/numSamples;
tic;
for i=1:numSamples
    feature(i,1) = getTotalCurvature(FV{i}.vertices, d);
end
time(2) = toc/numSamples;
tic;
for i=1:numSamples
    feature(i,2) = getRelativeLength(FV{i}.vertices, d);
end
time(3) = toc/numSamples;
time
%tic;
%feature(i,3) = getMirrorSurface(FV{i}.vertices, d);
%time(i,3) = toc;
%%



fig(1) = figure(1);

nbins = 50;
figure(1);
subplot(2,1,1);
maxCount = 10000;

hist(feature(:,1)',nbins);
axis([min(feature(:,1)) max(feature(:,1)) 0 maxCount]);
grid on;
title(['Features of ' int2str(numSamples) ' random samples']);
ylabel('Curvature');

subplot(2,1,2);
hist(feature(:,2)',nbins);
axis([min(feature(:,2)) max(feature(:,2)) 0 maxCount]);
grid on;
ylabel('Length');
%subplot(3,1,3);
%sel = feature(:,3) < 5e4;
%hist(feature(sel,3)',nbins);
%axis([0 5e4 0 300]);
%grid on;
%ylabel('MirrorSurface');


%save_figures(fig, './', ['feature_sampling_'], 12, [7 5]);
