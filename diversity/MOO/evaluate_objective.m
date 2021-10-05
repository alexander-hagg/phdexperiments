function fitness = evaluate_objective(genomes,d,varargin)
%OBJECTIVEMOOVSQD All objectives and features in one function for NSGA-II
if nargin == 0; fitness = 3; return; end % Return number of objective functions
invert = false; if nargin > 2; invert = varargin{1}; end

[fitness,polygons,~] = d.fitfun(genomes);
fitness = fitness';
features = categorize(polygons, d)';

% Area and Perimeter
if numel(polygons) > 1
    if invert
        fitness(:,2) = features(1,:); % max area
        fitness(:,3) = features(2,:); % max perimeter
    else
        fitness(:,2) = features(1,:); % max area
        fitness(:,3) = 1-features(2,:); % min perimeter
    end
else
    if invert
        fitness(2) = features(1); % max area
        fitness(3) = features(2); % max perimeter
    else
        fitness(2) = features(1); % max area
        fitness(3) = 1-features(2); % min perimeter
    end
end

end

