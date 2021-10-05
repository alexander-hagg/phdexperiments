function penalty = constraintPenalty(X,constraints)
%CONSTRAINTPENALTY Calculate constraint penalty
% penalty between 0 and 1
penalty = zeros(size(X,1),1);

[~,classDistances] = applyConstraints(X, constraints);
classBinary = logical(zeros(size(classDistances,1),1));
classBinary(constraints.selectedClasses) = 1;
distSELECTED = min(classDistances(classBinary,:));
distNONSELECTED = min(classDistances(~classBinary,:));

if ~isempty(distNONSELECTED) % If not all classes are selected    
    penalty = distSELECTED./(distSELECTED+distNONSELECTED); 
    penalty = reshape(penalty,length(penalty),1); % Make sure it is a col vector    
end
end

