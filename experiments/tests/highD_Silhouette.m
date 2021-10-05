
dims = 1:20;
numSamples = 5000;

for i=1:length(dims) 
    mu = zeros(1,dims(i));mu(1)=1;
    sigma = eye(dims(i));
    R = chol(sigma);
    z1 = zeros(numSamples,dims(i)) + randn(numSamples,dims(i))*R;
    z2 = repmat(mu,numSamples,1) + randn(numSamples,dims(i))*R;
    locations = [z1;z2];
    clusterLabels = [zeros(numSamples,1);ones(numSamples,1)];
    val(i) = m_silhouette(locations,clusterLabels);
end
%%
fig(1) = figure(1);
p = plot(val);
grid on;
xlabel('Dimensionality');
ylabel('Silhouette Index');
axis([0 max(dims) 0 0.25]);
p.LineWidth = 2;
save_figures(fig, './', ['silhouette_1-20D'], 20, [6 5]);

