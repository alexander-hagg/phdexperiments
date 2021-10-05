x = randn(1000,dims);

no_dims = 2; initial_dims = dims; theta = 0.5; alg = 'svd'; max_iter = 1000;
perplexity1= 5;
reducedParSpace1 = fast_tsne(x, no_dims, initial_dims, perplexity1, theta, alg, max_iter);


clear fig;

fig(1) = figure(1);hold off;
h1 = scatter(x(:,1)', x(:,2)','filled');
title(['Gaussian (orig. 1st 2 dims)']);
grid on;axis equal;
drawnow;

fig(2) = figure(2);hold off;
h1 = scatter(reducedParSpace1(:,1)', reducedParSpace1(:,2)','filled');
title(['Gaussian (reduced dims), perplexity ' num2str(perplexity1)]);
grid on;axis equal;
drawnow;


no_dims = 2; initial_dims = dims; theta = 0.5; alg = 'svd'; max_iter = 1000;
perplexity2= 50;
reducedParSpace2 = fast_tsne(x, no_dims, initial_dims, perplexity2, theta, alg, max_iter);

fig(3) = figure(3);hold off;
h1 = scatter(reducedParSpace2(:,1)', reducedParSpace2(:,2)','filled');
title(['Gaussian (reduced dims), perplexity ' num2str(perplexity2)]);
grid on;axis equal;
drawnow;
save_figures(fig, './', ['tsne_perplexity_'], 12, [7 7]);


