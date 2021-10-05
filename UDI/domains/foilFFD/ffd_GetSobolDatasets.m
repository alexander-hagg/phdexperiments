function filenames = ffd_GetSobolDatasets(d, varargin)
parse = inputParser;
parse.addRequired('d');
parse.addOptional('rndSampleSizes', [500, 1000]);
parse.parse(d, varargin{:});
rndSampleSizes = parse.Results.rndSampleSizes;

%d.express, d.categorize, d.featureMin, d.featureMax
getCoordinates = @(x) feval(d.categorize, x, d);
edges = {0:0.04:1 0:0.04:1};

filenames = string(rndSampleSizes);
filenames = strcat(filenames,'.mat');
filenames = strcat('ffd_',filenames);

for r=1:length(rndSampleSizes)
    sobSequence             = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
    
    observation             = nan(rndSampleSizes(r), d.dof); % new values will be stored here
    noval = any(isnan(observation'))';
    sobPoint = 1;
    while(sum(noval)>0)
        disp(['Total NaNs now: ' int2str(sum(noval))]);
        % Get observation
        observation(noval,:) = sobSequence(sobPoint:-1+sobPoint + sum(noval),:);
        % Get Map Indices
        mapLinIndx(noval)    = getBins( observation(noval,:), getCoordinates, d, edges);
        % Get Fitness Values
        [trueVal, trueFit]   = feval(d.preciseEvaluate, observation(noval,:), d);
        cD_true(noval)       = trueVal(:,1);
        cL_true(noval)       = trueVal(:,2);
        fitness_true(noval)  = trueFit';
        noval = isnan(cD_true');
        sobPoint = sobPoint + sum(noval);
    end
    save(filenames{r}, 'observation', 'mapLinIndx', 'cD_true', 'cL_true', 'fitness_true', '-v7.3');
end
end