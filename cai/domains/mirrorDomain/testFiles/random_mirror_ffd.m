d.dof = 36;

mutation = 0.5*rand(1,d.dof);
[FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');
view(90,0);
axis([-100 100 -200 200 600 800]);    

%% Create mirrors from Sobol set
p.numInitSamples = 100;
sobSequence                         = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
samples                             = sobSequence(1:(1+p.numInitSamples)-1,:);


for i=1:p.numInitSamples
    %[FV, ~, ffdP] = mirror_ffd_Express(0.5*ones(1,d.dof), 'mirrorBase.stl');
    [FV, ~, ffdP] = mirror_ffd_Express(samples(i,:), 'mirrorBase.stl');
    stlwrite(['test_' int2str(i) '.stl'],FV);

%    view(90,0);
%    axis([-100 100 -200 200 600 800]);   
%    pause(0.5);
end