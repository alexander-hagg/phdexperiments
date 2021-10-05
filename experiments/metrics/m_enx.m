function val = m_enx(distances_orignal,distances_reduced)
%ENX Summary of this function goes here
%   Detailed explanation goes here
k = 5; 
n = size(distances_orignal,1);
K_threshold = 15;
K_cutoff = 400;

coranking = zeros(size(distances_orignal));
[~,orderDd] = sort(distances_orignal);
[~,orderDq] = sort(distances_reduced);
corankentries = cat(3,orderDd,orderDq);
for ii=1:size(corankentries,1)
    for jj=1:size(corankentries,2)
        coranking(corankentries(ii,jj,1),corankentries(ii,jj,2)) = coranking(corankentries(ii,jj,1),corankentries(ii,jj,2)) + 1;
    end
end
% fig(i) = figure(i);imagesc(coranking);xlabel('low-dimensional space');ylabel('high-dimensional space');cc = colorbar;colormap(flipud(gray));cc.Label.String = '# Occurrences';
val = 0;
for kk=-K_threshold:K_threshold
    val = val + sum(sum(diag(coranking,kk)));
end
val = 1/(K_cutoff*n)*val;
end
