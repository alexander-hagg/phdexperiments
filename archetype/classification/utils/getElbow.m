function [maxVal, maxID, maxDist] = getElbow( yvals )
%GETELBOW Summary of this function goes here
%   Detailed explanation goes here

    % get coordinates of all the points
    nPoints = length(yvals);
    allCoord = [[1:nPoints]',yvals];              %'# SO formatting    
    firstPoint = allCoord(1,:);    
    % get vector between first and last point - this is the line
    lineVec = allCoord(end,:) - firstPoint;    
    % normalize the line vector
    lineVecN = lineVec / sqrt(sum(lineVec.^2));    
    % find the distance from each point to the line:
    vecFromFirst = bsxfun(@minus, allCoord, firstPoint);
    
    % To calculate the distance to the line, we split vecFromFirst into two
    % components, one that is parallel to the line and one that is perpendicular
    % Then, we take the norm of the part that is perpendicular to the line and
    % get the distance.
    % We find the vector parallel to the line by projecting vecFromFirst onto
    % the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
    % We project vecFromFirst by taking the scalar product of the vector with
    % the unit vector that points in the direction of the line (this gives us
    % the length of the projection of vecFromFirst onto the line). If we
    % multiply the scalar product by the unit vector, we have vecFromFirstParallel
    scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
    vecFromFirstParallel = scalarProduct * lineVecN;
    vecToLine = vecFromFirst - vecFromFirstParallel;
    
    % Make sure we do not take an upper elbow.
    vecToLine(vecToLine>0) = 0;
    
    % distance to line is the norm of vecToLine
    distToLine = sqrt(sum(vecToLine.^2,2));
    [maxDist, maxID] = max(distToLine);
    maxVal = yvals(maxID);
    
end

