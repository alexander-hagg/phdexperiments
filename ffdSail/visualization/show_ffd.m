clear fig;
fig(1) = figure(99);hold off;
base = loadBaseAirfoil(allresults{1,1}.d.express, 10);
plot(base.foil(1,:),base.foil(2,:), 'Color', [0 0 0], 'LineWidth',2);
hold on;
drawnow;
genome= 0.5*ones(1,10);
for g=1:5
    foil = ffdRaeY(genome);
    plot(foil(1,:),foil(2,:), 'Color', [g/5 0 0], 'LineWidth',2);
    axis equal;
    axis tight;
    drawnow;
    pause(0.01);
    genome(2) = genome(2)+0.2;    
end
l = legend({'Base Foil'});
l.FontSize = 14;
save_figures(fig, './', ['ffd_'], 14, [12 12]);
