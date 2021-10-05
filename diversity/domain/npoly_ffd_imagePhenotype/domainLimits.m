% Find feature bounds for MAP-Elites
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = p.numInitSamples+1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+100000)-1,:)+d.ranges(:,1)';

[fitness,polygons] = d.fitfun(initSamples);
[~,nonNormalizedFeatures] = categorize(polygons, d);
measuredFeatureMins = min(squeeze(nonNormalizedFeatures));
measuredFeatureMaxs = max(squeeze(nonNormalizedFeatures));

%%
measuredFeatureMins
measuredFeatureMaxs
figure(1);
for i=1:length(polygons)
    plot(polygons{i});
    title([mat2str(nonNormalizedFeatures(i,:)) ' - ' int2str(i)]);
    drawnow;
    pause(2)
end
%%
figure(1); 
    
for i=1:size(initSamples,1)
    subplot(5,5,i);
    plot(polygons{i})
    fitnessPointSymmetry(initSamples(i,:),d)
    axis([-0.5 0.5 -0.5 0.5]);
    aa = area(polygons{i});
    pp = perimeter(polygons{i});
    title([num2str(aa) ' - ' num2str(pp)]);
    drawnow;
    
    pause(0.5)
end

%%
circle = [ones(1,DOF/2),pi*ones(1,DOF/2)]
[fit,pp] = fitnessPointSymmetry(circle,d);
figure(1);
plot(pp{1}); axis([-0.5 0.5 -0.5 0.5]);
area(pp{1})
perimeter(pp{1})

star = [[1 0.01 1 0.01 1 0.01 1 0.01],0*ones(1,DOF/2)]
[fit,pp] = fitnessPointSymmetry(star,d);
figure(2);
plot(pp{1}); axis([-0.5 0.5 -0.5 0.5]);
area(pp{1})
perimeter(pp{1})

