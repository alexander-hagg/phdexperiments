function predMap = getGroundTruthMapArrayMirror(predMap, d, varargin)
%% Get true results of each prediction map for the mirror case
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
        trueVal = [predMap(i).dragForce(:)]./10;
    else
        [trueVal, trueFit] = feval(d.preciseEvaluate,genes,d);
    end
    
    % Initialize the maps in the correct shape
    predMap(i).fitness_true         = nan(size(predMap(i).fitness));
    predMap(i).dragForce_true       = nan(size(predMap(i).dragForce));
    
    % Insert Values
    predMap(i).fitness_true(:)      = trueFit;
    predMap(i).dragForce_true(:)    = trueVal(:,1);
    
    disp(num2str(toc));
end

disp(['Run done.'])

%------------- END OF CODE --------------






