%% More Components is better


figs(1) = figure(1);
hold off;
X = randn(50,2);
scatter(X(:,1),X(:,2),32,[0 0 0],'filled');
hold on;
X = randn(50,2);
scatter(X(:,1),0.05*X(:,2),32,[1 0 0],'filled');

grid on;
boundsX = 2;
boundsY = 2;
axis([-boundsX boundsX -boundsY boundsY]); 
title('Principal Components');
legend('More Diverse','Less Diverse');
xlabel('PC1');
ylabel('PC2');

%% Larger is better
X = randn(10,2);

scaleY = 0.5;

figs(2) = figure(2);
hold off;
scatter(X(:,1),scaleY*X(:,2),32,[0 0 0],'filled');
hold on;
scatter(0.2*X(:,1),0.2*scaleY*X(:,2),32,[1 0 0],'filled');

plot([X(:,1) 0.2*X(:,1)]',[scaleY*X(:,2) 0.2*scaleY*X(:,2)]','Color','k');

grid on;
boundsX = 2;
boundsY = 2;
axis([-boundsX boundsX -boundsY boundsY]); 
title('Principal Components');
legend('More Diverse','Less Diverse');
xlabel('PC1');
ylabel('PC2');

%% Better is better
X = randn(50,3);
X(:,3) = X(:,3) - min(X(:,3));
figs(3) = figure(3);
hold off;
scatter(X(:,1),0.2+0.1*X(:,3),32,[0 0 0],'filled');
hold on;
scatter(X(:,1),0.1*X(:,3),32,[1 0 0],'filled');

plot([X(:,1) X(:,1)]', [0.2+0.1*X(:,3) 0.1*X(:,3)]','Color','k');
grid on;
axis([-2 2 0 1]); 
title('Principal Components');
legend('Higher Quality','Lower Quality');
xlabel('PC1');
ylabel('PC2');
zlabel('Fitness');
%%
save_figures(figs,'.','compressibility_diversity',14,[6 4])
