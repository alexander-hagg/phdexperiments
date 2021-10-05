ranges = d.ranges;
radialBoundAdaptation = 0.25;
axialBoundAdaptation = 0.01;
numInitSamples = 100000;
ranges(:,1) = [axialBoundAdaptation*ones(d.dof/2,1);-radialBoundAdaptation*pi*ones(d.dof/2,1)];
ranges(:,2) = [ ones(d.dof/2,1); radialBoundAdaptation*pi*ones(d.dof/2,1)];


sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen'); sobPoint = 1;
sobolSample = sobSequence(sobPoint:(sobPoint+numInitSamples)-1,:);
samples = range(ranges').*sobolSample+ranges(:,1)';

samples(end+1,:) = max(samples);
    %%
showPhenotype(samples,d)
%%
pgons = getPhenotypeFFD(samples,d.base);
for i=1:length(pgons)
    areas(i) = area(pgons{i});
    circumference(i) = perimeter(pgons{i});
end

disp(['Min Area: ' num2str(min(areas))]);
disp(['Max Area: ' num2str(max(areas))]);
disp(['Min circumference: ' num2str(min(circumference))]);
disp(['Max circumference: ' num2str(max(circumference))]);

%% Basic shapes
%[features,unnormalizedFeatures] = categorize(pgons,d)
%max(unnormalizedFeatures)
%ans =
%    0.7067    3.3997