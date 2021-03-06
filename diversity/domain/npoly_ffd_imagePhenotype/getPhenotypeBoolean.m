function [flatbooleanMap,booleanMap] = getPhenotypeBoolean(polyshapes,varargin)
%GETPHENOTYPEBOOLEAN Summary of this function goes here
%   Detailed explanation goes here

resolution = 64; if nargin>1;resolution = varargin{1};end

if isempty(polyshapes)
    flatbooleanMap = [];
    booleanMap = [];
    return;
end

for i=1:length(polyshapes)
        pixelCoordinates = ceil(polyshapes{i}.Vertices*(resolution)/2)+resolution/2;
        pixelCoordinates(all(isnan(pixelCoordinates)'),:) = [];
        booleanMap{i} = poly2mask(pixelCoordinates(:,1),pixelCoordinates(:,2),resolution,resolution);
end


flatbooleanMap = cat(3,[],booleanMap{:});
flatbooleanMap = reshape(flatbooleanMap,[],size(flatbooleanMap,3));
flatbooleanMap = flatbooleanMap';

end

