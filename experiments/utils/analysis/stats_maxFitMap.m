function maxFitMap = stats_maxFitMap(mapSequences)
%%MAXFITMAP Get maximally fit map from multiple sequences of feature maps
for j=1:length(mapSequences)
    t_maxFitMap = zeros(size(mapSequences{j}.predMap(end).fitness_true));
    for ii=1:length(mapSequences{j}.predMap)
        thisMap = mapSequences{j}.predMap(ii).fitness_true;
        betterValues = logical( (t_maxFitMap > thisMap) .* (~isnan(thisMap)) );
        t_maxFitMap(betterValues) = thisMap(betterValues);
    end
    if j==1
        maxFitMap = t_maxFitMap;
    else
        improvedBins = maxFitMap > t_maxFitMap;
        maxFitMap(improvedBins) = t_maxFitMap(improvedBins);
    end
end

end


