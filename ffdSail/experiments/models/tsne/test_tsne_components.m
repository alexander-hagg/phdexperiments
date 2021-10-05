set(0,'DefaultFigureWindowStyle','docked')

xsamples = [-2 0 8];

%%
fig(1) = figure(1);hold off;
scatter(xsamples(1),zeros(1,1), 64, [0 0 1], 'filled');
hold on;
scatter(xsamples(2),zeros(1,1), 64, [0 1 0], 'filled');
scatter(xsamples(3),zeros(1,1), 64, [1 0 0], 'filled');
axis([-20 20 -0.05 0.2]);
grid on;
ax = gca;
ax.XTick = xsamples;
ax.XTickLabel = strcat('x_',{'j','i','k'});
legend('Datapoint', 'Compare with', 'Other datapoint');
xlabel('High-Dimensional Input');
ylabel('pdf');
%title('Similarity comparison');
%%
fig(2) = figure(2);hold off;
x = [-20:.1:20];
plot(x,normpdf(x,0,4),'b','LineWidth',2);
hold on;

scatter(xsamples(1),zeros(1,1), 64, [0 0 1], 'filled');
scatter(xsamples(2),zeros(1,1), 64, [0 1 0], 'filled');
scatter(xsamples(3),zeros(1,1), 64, [1 0 0], 'filled');
plot([xsamples(2) xsamples(2)],[0 0.2], 64, [0 1 0], 'LineWidth',2);

axis([-20 20 -0.05 0.2]);
grid on;
ax = gca;
ax.XTick = xsamples;
ax.XTickLabel = strcat('x_',{'j','i','k'});
legend('Distr. around x_1', 'Datapoint', 'Compare with', 'Other datapoint');
xlabel('High-Dimensional Input');
ylabel('pdf');
%title('Distribution around x_i');

%%
fig(3) = figure(3);hold off;
x = [-20:.1:20];
plot(x,normpdf(x,0,4),'b','LineWidth',2);
hold on;

scatter(xsamples(1),zeros(1,1), 64, [0 0 1], 'filled');
scatter(xsamples(2),zeros(1,1), 64, [0 1 0], 'filled');
scatter(xsamples(3),zeros(1,1), 64, [1 0 0], 'filled');
plot([xsamples(1) xsamples(1)],[0 normpdf(xsamples(1),0,4)], 64, [0 0 1]);
plot([xsamples(2) xsamples(2)],[0 normpdf(xsamples(2),0,4)], 64, [0 1 0]);
plot([xsamples(3) xsamples(3)],[0 normpdf(xsamples(3),0,4)], 64, [1 0 0]);

axis([-20 20 -0.05 0.2]);
grid on;
ax = gca;
ax.XTick = xsamples;
ax.XTickLabel = strcat('x_',{'j','i','k'});
ax.YTick = [normpdf(xsamples(3),0,4) normpdf(xsamples(1),0,4) normpdf(xsamples(2),0,4)]
legend('Distr. around x_1', 'Datapoint', 'Compare with', 'Other datapoint');
xlabel('High-Dimensional Input');
ylabel('pdf');
%title('Conditional Probability');

%%
fig(4) = figure(4);hold off;
x = [-20:.1:20];
plot(x,normpdf(x,0,1),'b','LineWidth',2);
hold on;
plot(x,tpdf(x,1),'r','LineWidth',2);
grid on;
legend('Normal', 'Student-t');

%save_figures(fig, './', ['tsne_lecture_'], 12, [7 3]);

