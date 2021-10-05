function h = viewCellRanking( pVals, tVals, cat, edges, lstyle, clrs, lw, selection )
%VIEWCELLRANKING Summary of this function goes here
%   Detailed explanation goes here
for sampleSet=1:length(pVals)
    for modelSelIndex=1:length(selection)
        modelIndex=selection(modelSelIndex);
        for rep=1:size(pVals{sampleSet},2)
            validRankings(modelSelIndex, sampleSet, rep) = 0;
            allRankings(modelSelIndex, sampleSet, rep) = 0;
            allPredictions  = pVals{sampleSet}(modelIndex,rep,:,1);
            allTrues        = tVals{sampleSet}(modelIndex,rep,:,:);
            % Assign feature bins
            feature         = squeeze(cat{sampleSet}(rep,:,:));
            clear bin; for iDim = 1:length(edges); bin(:,iDim) = discretize(feature(:,iDim),edges{iDim}); end
            [~,~,ic] = unique(bin(:,[1 2]),'rows');
            
            % All comparisons within bins
            for i=1:max(ic)
                t_predictions = squeeze(allPredictions(ic==i));
                if numel(t_predictions) > 1
                    cmpmat = repmat(t_predictions,1,length(t_predictions));
                    PRANK = cmpmat < t_predictions';
                    
                    t_trues = squeeze(allTrues(ic==i));
                    cmpmat = repmat(t_trues,1,length(t_trues));
                    TRANK = cmpmat < t_trues';
                    validRankings(modelSelIndex, sampleSet, rep) = validRankings(modelSelIndex, sampleSet, rep) + sum(PRANK(:)&TRANK(:)) + sum(~PRANK(:)&~TRANK(:));
                    allRankings(modelSelIndex, sampleSet, rep) = allRankings(modelSelIndex, sampleSet, rep) + numel(PRANK);
                end
            end
        end
    end
end
%%
relRanks = 100 - 100*validRankings./allRankings;
for modelIndex=1:size(relRanks,1)
    %le = prctile(relRanks(modelIndex,:,:),25,3);le(1) = [];
    %ri = prctile(relRanks(modelIndex,:,:),75,3);ri(1) = [];
    %patch([2:length(le)+1,length(le)+1:-1:2],[ri fliplr(le)], clrs(modelIndex,:),'FaceAlpha',.5, 'EdgeColor', clrs(modelIndex,:), 'LineWidth', 2);
    h(modelIndex) = plot(median(relRanks(modelIndex,:,:),3), lstyle{mod(modelIndex,2)+1}, 'LineWidth', lw, 'Color', clrs(modelIndex,:));
    hold on;
    
end
%%
end

