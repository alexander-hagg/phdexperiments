function predMap = getGroundTruthMapArray(predMap, d, varargin)
%% Get true results of each prediction map
%

%------------- BEGIN CODE --------------
% Read in input
parse = inputParser;
parse.addRequired('map');
parse.addRequired('d');
%parse.addOptional('numSamples',[]);
parse.addOptional('dummy',false);
parse.parse(predMap, d, varargin{:});
%iItr = parse.Results.numSamples; 
dummy = parse.Results.dummy;

tic;
for i = 1:length({predMap.fitness})
    if isempty(predMap(i).fitness); continue;end
    
    disp(['Iteration ' int2str(i) '/' int2str(length({predMap.fitness}))]);
    genes = reshape(predMap(i).genes,...
        [size(predMap(i).genes,1) * ...
        size(predMap(i).genes,2), d.dof]) ;
    
    % Just test if it's working?
    if dummy
        trueFit = predMap(i).fitness(:)./10;
        trueVal = [predMap(i).cD(:) predMap(i).cL(:)]./10;
    else
        [trueVal, trueFit] = feval(d.preciseEvaluate,genes(~any(isnan(genes')),:),d);
    end
    
    % Initialize the maps in the correct shape
    predMap(i).fitness_true    = nan(size(predMap(i).fitness));
    predMap(i).cD_true         = nan(size(predMap(i).cD));
    predMap(i).cL_true         = nan(size(predMap(i).cL));
    
    % Insert Values
    predMap(i).fitness_true(~any(isnan(genes'))) = trueFit;
    predMap(i).cD_true(~any(isnan(genes')))      = trueVal(:,1);
    predMap(i).cL_true(~any(isnan(genes')))      = trueVal(:,2);
    
    disp(num2str(toc));
end

disp(['Run done.'])

%------------- END OF CODE --------------






