function VIS_allrobots(runData,itID)
%VIS_ALLROBOTS Summary of this function goes here
%   Detailed explanation goes here
%%
clear isOutside
genes = reshape(runData{itID}.acqMap.genes,900,size(runData{itID}.acqMap.genes,3));
genes = genes(all(~isnan(genes')),:);
d = runData{itID}.d;
p = runData{itID}.p;
clear goalID dist
%trajectories = runData{itID}.acqMap.phenotype;
d.tmpdir = '/tmp';
[trajectories] = eval_maze(genes,d);
for i=1:size(trajectories,1)
    [goalID(i),dist(i,:,:)] = evalFirstGoal(squeeze(trajectories(i,:,:)),d.goalLocations);
end
clrs = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
%
%isOutside = any(squeeze(dist(:,4,:))' > 336/2); trajectories = trajectories(isOutside,:,:);
if exist('isOutside'); goalID = goalID(isOutside);end
disp('done analysis');
%%
writeFig = true;
filename = ['robotsSelOUT3_theta_' int2str(runData{itID}.d.theta) '_it_' int2str(itID) '_selVal_' num2str(runData{itID}.p.selectionValue) '.gif'];
delaytime = 0.1;startFrame = 10;frameJump = 100;
showPath = true;
goalNames = {'red','green','blue','none'};
selValue = runData{itID}.p.selectionValue; if selValue > 10; selValue = 4; end
h = figure;
for endDraw = [startFrame:frameJump:size(trajectories,2) size(trajectories,2)]
    %%
    hold off;
    
    drawMap = true;plotSim(squeeze(trajectories(1,1:endDraw,:)),d,drawMap,[0 0 0],showPath);
    drawMap = false;plotSim(trajectories(:,1:endDraw,:),d,drawMap,clrs(goalID+1,:),showPath);
    drawnow;
    title(['Exit selected: ' goalNames{selValue} ', Iteration: ' int2str(itID) ', Time Step: ' int2str(endDraw)]);
    if writeFig
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if endDraw == startFrame
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delaytime);
        elseif (size(trajectories,2)-endDraw)< frameJump
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',10*delaytime);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delaytime);
        end
    end
end

end

