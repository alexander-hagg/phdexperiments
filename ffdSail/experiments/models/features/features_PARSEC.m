%% Get samples
[samples, value] = initialSampling(d,2000);


%% Get Features and Visualize
shape = d.express(samples);
[xup,idxup] = max(squeeze(shape(2,1:end/2,:)));
xup = xup'; idxup = idxup';
idxup = idxup/200;
curvature = squeeze(sum(abs(shape(2,2:end,:)-shape(2,1:end-1,:))));
maxthickness = squeeze(max(abs(shape(2,1:end/2,:)-fliplr(shape(2,end/2+1:end,:)))));
area = squeeze(polyarea(shape(1,:,:), shape(2,:,:)));
camber = squeeze(shape(2,1:end/2,:)+shape(2,end:-1:end/2+1,:))/2;
maxcamber = max(camber)';
%%
set(0,'DefaultFigureWindowStyle','default')
fig(1) = figure(1);
plot(area,curvature,'x');xlabel('total area');ylabel('total curvature');
grid on
axis equal;
axis([0 1 0 1]);
corr = corrcoef(area,curvature);
title(['Corr. Coeff. ' num2str(corr(2))])

fig(2) = figure(2);
plot(maxthickness,curvature,'x');xlabel('max thickness');ylabel('total curvature');
grid on
axis equal;
axis([0 1 0 1]);
corr = corrcoef(maxthickness,curvature);
title(['Corr. Coeff. ' num2str(corr(2))])

fig(3) = figure(3);
plot(xup,curvature,'x');xlabel('xup');ylabel('total curvature');
grid on
axis equal;
axis([0 1 0 1]);
corr = corrcoef(xup,curvature);
title(['Corr. Coeff. ' num2str(corr(2))])

fig(4) = figure(4);
plot(idxup,curvature,'x');xlabel('location xup');ylabel('total curvature');
grid on
axis equal;
axis([0 1 0 1]);
corr = corrcoef(idxup,curvature);
title(['Corr. Coeff. ' num2str(corr(2))])

fig(5) = figure(5);
plot(maxcamber,curvature,'x');xlabel('max camber');ylabel('total curvature');
grid on
axis equal;
axis([0 1 0 1]);
corr = corrcoef(maxcamber,curvature);
title(['Corr. Coeff. ' num2str(corr(2))])

save_figures(fig, './', ['features_correlation_'], 12, [5 5]);

%% Correlation features/parameters

features = [area, curvature, maxthickness, xup, idxup, maxcamber];

for f=1:size(features,2)
    for p=1:size(samples,2)
        cr = corrcoef(features(:,f),samples(:,p));
        corrFeatPar(f,p) = cr(2);
    end
end
max(abs(corrFeatPar'))
mean(abs(corrFeatPar'))

for f=1:size(features,2)
    for p=1:size(features,2)
        cr = corrcoef(features(:,f),features(:,p));
        corrFeatFeat(f,p) = cr(2);
    end
end
max(corrFeatFeat')

%% Picked features

fig(1) = figure(1);
plot(area,maxcamber,'x');xlabel('total area');ylabel('max camber');
grid on
axis equal;
axis([0 0.2 0 0.2]);
corr = corrcoef(area,maxcamber);
title(['Corr. Coeff. ' num2str(corr(2))])

save_figures(fig, './', ['features_correlation_'], 12, [5 5]);
