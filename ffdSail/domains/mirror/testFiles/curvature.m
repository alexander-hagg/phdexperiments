[FV, ~, ffdP] = mirror_ffd_Express(0.5 + zeros(1,41), 'mirrorBase.stl');

% Visualize
fig(1) = figure(2);hold off;
plot3(FV.vertices(1,:), FV.vertices(2,:), FV.vertices(3,:),'x');
hold on;
plot3(FV.vertices(1,baseSubMeshIds), FV.vertices(2,baseSubMeshIds), FV.vertices(3,baseSubMeshIds),'x');
view(112,90);

lineIDs{1} = [8451 8107 7633 7109 6421 5603 5109 4761 4433 4090 3788 3568 3378 3204 3064 2920 2788 2676 2582 2514 2460 2464 2495 2524 2620 2792 2962 3127 3305 3476 3624 3910 4172 4419 4659 4870 5082];
lineIDs{2} = [8433 8134 7771 7347 6831 6037 5437 5024 4614 4222 3943 3711 3527 3279 3100 2949 2835 2734 2652 2573 2508 2420 2335 2259 2204 2174];
lineIDs{3} = [8022 7755 7400 6970 6446 5835 5437 5205 4952 4717 4761 4945 5214 5421 5689 6017 6413 6997 7667];
lineIDs{4} = [7544 6942 6579 5292 4749 4283 3979 3788 3653 3570 3517 3527 3486 3525 3619 3776 3998 4322 4604 4664 4730 4826 4917];
lineIDs{5} = [7232 5535 4862 4077 3527 3310 3185 3111 3064 3050 3092 3164 3309 3532 3893 4151 4375 4584 4797 5007 5211];
lineIDs{6} = [8526 8507 8504 8531 8433 8354 8353 8292 8246 8022 8005 7911 7816 7737 7661 7605 7385 7429 7321 7232 6875 6872 6762 6232 6022 5953 6001 5738 5554 5456 5329 5290 5063 5299 5126 4917];
lineIDs{7} = [8507 8526 8472 8434 8461 8506 8459 8451 8510 8413 8329 8274 8226 8186 8190 8014 8113 8039 7845 7896 7676 7559 7435 7326 7245 7065 6944 6866 6678 6605 6437 6281 6178 5797 5793 5807 5605 5696 5564 5447 5353 5272 5211 5159 5113 5082];


for i=1:length(lineIDs)
    pp = FV.vertices(:,lineIDs{i});
    plot3(pp(1,:),pp(2,:),pp(3,:),'ko','LineWidth',8);
end

for i=1:length(lines)
    scatter3(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),'filled');
    for ii=1:length(lines{i})
        [~,ID] = intersect(FV.vertices',lines{i}(ii,:),'rows');
        text(lines{i}(ii,1),lines{i}(ii,2)+0.7,lines{i}(ii,3), int2str(ID));
    end
end
%

xlabel('x');ylabel('y');zlabel('z');
%view(-270,0);
grid on;

%% Find IDs
[~,index_InMesh] = intersect(FV.vertices',back,'rows')
for i=1:length(lines)
    [~,indices] = intersect(back,lines{i},'rows')
    linIndices{i} = index_InMesh(indices)'
end
indexLine = cell2mat(linIndices)
%% Visualize random shapes and curvature
% Select planes
%if ~exist('curvs','var')
    for ii=1:1000
        [FV, ~, ffdP] = ...
            mirror_ffd_Express(0.5 + 0.2*randn(1,41), 'mirrorBase.stl');
        % figure(3);
        % hold off;
        % plot3(FV.vertices(1,:),FV.vertices(2,:),FV.vertices(3,:),'x');
        % hold on;
        % for i=1:5
        % scatter3(FV.vertices(1,linIndices{i}),FV.vertices(2,linIndices{i}),FV.vertices(3,linIndices{i}),'filled');
        % end
        % %
        tic
        totalCurv(ii) = 0;
        for i=1:length(linIndices)
            line = FV.vertices(:,linIndices{i});
            % Total sum of absolute curvatures
            planes = logical([1 1 0; 1 0 1; 0 1 1]);
            for plane=1:size(planes,1)
                totalCurv(ii) = totalCurv(ii) + sum(abs(LineCurvature2D(line(planes(plane,:),:)')));
            end
        end
        toc
        for i=1:length(linIndices)
            %line = FV.vertices(:,linIndices{i});
            %     text(line(1,1),line(2,1),line(3,1),num2str(curv(i)),'FontSize',14)
        end
        title(['Total Curvature Approximation: ' num2str(totalCurv(ii))]);
        disp(['Total Curvature Approximation: ' num2str(totalCurv(ii))]);
       
    end
%end
figure(4);
histogram(totalCurv);xlabel('Total Curvature');ylabel('pdf');title('Curvature of 1000 random samples');grid on;

%% Get some random shapes and show curvature
fig(1) = figure(5);
for ii=1:9
    subplot(3,3,ii);
    hold off;
    [FV{ii}, ~, ffdP{ii}] = ...
        mirror_ffd_Express(0.5 + ii*0.05*randn(1,41), 'mirrorBase.stl');
    scatter3(FV{ii}.vertices(1,:), FV{ii}.vertices(2,:),FV{ii}.vertices(3,:),4,'filled');
    view(-30,10);
    hold on;
    tic
    for i=1:length(linIndices)
        line = FV{ii}.vertices(:,linIndices{i});
        % Ratio between line length and length between start and end point
        startToEnd(i) = pdist2(line(:,1)',line(:,end)');
        curv(i) = sum(diag(squareform(pdist(line')),1))/startToEnd(i);
        scatter3(line(1,:),line(2,:),line(3,:),'filled');
    end
    toc
    for i=1:length(linIndices)
        curv(i) = curv(i)/max(startToEnd);
    end
    title(['' num2str(sum(curv))]);
    disp(['Total Curvature Approximation: ' num2str(sum(curv))]);
    drawnow;
    
end
%%
save_figures(fig, './', ['fff_'], 12, [12 9]);