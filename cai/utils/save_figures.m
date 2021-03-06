function save_figures(fig, path, name, fontsize, canvas_size)
%SAVE_FIGURES Summary of this function goes here
% Save figures
for i = 1:length(fig)
    if isgraphics(fig(i))
        set(fig(i), 'Units', 'Centimeters');
        for c=1:length(fig(i).Children)
            if isa(fig(i).Children(c),'matlab.graphics.axis.Axes')
                fig(i).Children(c).FontSize = fontsize;
                fig(i).Children(c).XLabel.FontSize = fontsize;
                fig(i).Children(c).YLabel.FontSize = fontsize;
            end
        end
        fig(i).CurrentAxes.FontSize = fontsize;
        fig(i).CurrentAxes.XLabel.FontSize = fontsize;
        fig(i).CurrentAxes.YLabel.FontSize = fontsize;
        %savefig([path '/' name int2str(i) '.fig']);
        save2pdf([path '/' name int2str(i) '.pdf'],fig(i),600,canvas_size);
        system(['pdfcrop ' [path '/' name int2str(i) '.pdf'] ' ' [path '/' name int2str(i) '.pdf']]);
    end
end

end

