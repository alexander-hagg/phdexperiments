function fig = viewMirrors(d, optima, viewAngles, fontsize, rotation, figStartID, annotate, sizes)
%%VIEWMIRRORSMAP shows all mirror prototypes
%
%

if nargin < 5
    rotation = 90;
end
if nargin < 6
    figStartID = 1;
end
if nargin < 7
    annotate = false;
end

for i=1:size(optima,1)
    mutation = optima(i,:);
    %pred_fitness = prototypeFitness(i,j,:);
    [FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');
    
    for v=1:length(viewAngles)
        fig((i-1)*length(viewAngles)+v) = figure(figStartID + (i-1)*length(viewAngles)+v);
        delete(findall(gcf,'type','annotation'))
        hold off;
        hm = mirrorVisPaper(FV,ffdP, d, rotation, viewAngles(v), false, false);
        if hm~=0
            hm.EdgeColor = [0 0 0];
            hm.FaceColor = [1 1 1];
            hm.FaceAlpha = 1;
        end
        xlabel('x');ylabel('y');zlabel('z');
        view(rotation,viewAngles(v));grid on;axis equal;axis tight;
        ax = gca;
        %set(gca, 'visible', 'off')
        axis([-120 100 -150 120 -1000 1000]);
        
        

        %ax.Visible = 'off';
        if annotate
            str = ['Prototype ' int2str(i) ' with ' int2str(sizes(i)) ' members'];        
            dim = [.25 .5 .3 .3];
            a = annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
            a.FontSize = fontsize;
        end
        
    end
    
    %drawnow;
    
end

end