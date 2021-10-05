function val = m_gplus(distances,labels)
%GPLUS Summary of this function goes here
%   Detailed explanation goes here

numDiscordantPairs = 0;
numConcordantPairs = 0;
for ss=1:size(distances,1)
    label = labels(ss);
    thiscluster = find(labels==label);
    others = find(labels~=label);
    thisClusterDistances = distances(ss,thiscluster);
    for tt=1:length(thisClusterDistances)
        numDiscordantPairs = numDiscordantPairs + sum(thisClusterDistances(tt) > distances(ss,others));
        numConcordantPairs = numConcordantPairs + sum(thisClusterDistances(tt) <= distances(ss,others));
    end
end
t = (size(distances,1)*size(distances,1)-1)/2;
val = (2*numDiscordantPairs)/(t*(t-1));
end

