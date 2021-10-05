close all;set(0,'DefaultFigureWindowStyle','default');
%set(0,'DefaultFigureWindowStyle','docked')

clear figs;
scale = [7 6];
angles = [20];
rotation = 260;%80;

%prototypes = 0.5*ones(1,48);
%figs = viewMirrorPrototypes(prototypes, angles, 16, false);
%drawnow;%pause(1);
%save_figures(figs, './', ['PRODUQD_Run1_It0_'], 16, scale);

%
datastruct = output;filename = 'PRODUQD_SAIL';
%datastruct = prodigi_output{1}{1};filename = 'PRODUQD_Run1_It';
clear figs;

show = 1:5;
for i=1:length(datastruct)
    show(end+1) = datastruct{i}.conceptSelection.id;
    prototypes = datastruct{i}.prototypes(show,:);
    figs = viewMirrorPrototypes(prototypes, angles, 16, rotation, false);
    %drawnow;%pause(1);
    %save_figures(figs, './', ['PRODUQD_Run1_It' int2str(i) '_'], 16, scale);
end
