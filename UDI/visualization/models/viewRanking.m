function h = viewRanking( pVals, tVals, lstyle, clrs, lw, selection )
%VIEWRANKING Summary of this function goes here
%   Detailed explanation goes here
for sampleSet=1:length(pVals)
    for modelSelIndex=1:length(selection)
        modelIndex=selection(modelSelIndex);
        pr{sampleSet}(:,:,:) = squeeze(pVals{sampleSet}(modelIndex,:,:,1));
        cmpmat = repmat(pr{sampleSet}(:),1,length(pr{sampleSet}(:)));
        PRANK = cmpmat<pr{sampleSet}(:)';
        tr{sampleSet}(:,:,:) = squeeze(tVals{sampleSet}(modelIndex,:,:));
        cmpmat = repmat(tr{sampleSet}(:),1,length(tr{sampleSet}(:)));
        TRANK = cmpmat<tr{sampleSet}(:)';
        validRankings = sum(PRANK(:)&TRANK(:)) + sum(~PRANK(:)&~TRANK(:));
        relValidRankings(modelSelIndex,sampleSet) = 100 - 100*validRankings'/(size(PRANK,1)*size(PRANK,2));
    end
end

for modelIndex=1:size(relValidRankings,1)
    h(modelIndex) = plot(relValidRankings(modelIndex,:)', lstyle{mod(modelIndex,2)+1},'LineWidth', lw, 'Color', clrs(modelIndex,:));
    hold on;
end

end

