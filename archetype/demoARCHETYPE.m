% demoARCHETYPE - rapid architecture prototyping (without surrogate
% assistance)
%
% Computer Aided Ideation: Prototype Discovery using Quality-Diversity
%
% Please include the following references in any publication using this code. For Bibtex please see the end of this file.
%
%Hagg, A., Asteroth, A. and B?ck, T., 2018, September. Prototype discovery using quality-diversity. In International Conference on Parallel Problem Solving from Nature (pp. 500-511). Springer, Cham.
%Hagg, A., Asteroth, A. and B?ck, T., 2019, July. Modeling user selection in quality diversity. In Proceedings of the Genetic and Evolutionary Computation Conference (pp. 116-124). ACM.
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Oct 2019; Last revision: 01-Oct-2019

%------------- BEGIN CODE --------------
%% Configuration
clear;clc;
DOF = 16;
DOMAIN = 'npoly_ffd_archetype';
ALGORITHM = 'grid';
nIters = 2;
nPrototypes = 20;

addpath(genpath('.'));
iter = 1; app.selectedPrototypes = {}; app.constraints = {}; app.d = {}; app.p = {};
rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
app.d{iter} = domain(DOF);

rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);
app.p{iter} = defaultParamSet(4);

%% Main loop
for iter=1:nIters
    disp(['Iteration: ' int2str(iter)]);
    
    %% I) Illumination with QD
    disp(['Illumination']);
    [map,fitnessFunction] = initialize(app.constraints,app.d{iter},app.p{iter},false);
    app.map{iter} = illuminate(map,fitnessFunction,app.p{iter},app.d{iter});

                
    %% II) Extract prototypes
    disp(['Similarity space and prototype extraction']);
    [app.classification{iter}, simspace] = extractClasses(app.map{iter}.genes,[],'kmedoids',nPrototypes);
    % Show classes in similarity space
    f=figure(4);clf(f);figHandle=axes;cla(figHandle);viewClasses(app.map{iter}.genes,app.classification{iter},figHandle);title(figHandle,'Similarity Space');
    % Show prototypes
    numRows = ceil( size(app.classification{iter}.protoX,1).^(2/3) );
    for i=1:size(app.classification{iter}.protoX,1)
        placement(i,:) = [mod(i-1,numRows) floor((i-1)/numRows)]*app.d{iter}.spacer;
    end
    f=figure(5);clf(f);figHandle=axes;cla(figHandle);showPhenotype(figHandle,app.d{iter},app.classification{iter}.protoX,placement);title(figHandle,'Prototypes');
    
    %% III) Selection procedure
    prompt = 'Please select a prototype by entering an integer ID (bottom left to upper right: ';
    app.selectedPrototypes{iter} = input(prompt);
    
    %% IV) Create UDHM based on selected prototypes
    if length(app.selectedPrototypes) >= iter && ~isempty(app.selectedPrototypes{iter})
        app.constraints{iter} = setConstraints(app.classification{iter},app.selectedPrototypes{iter},[],'class');
        % Initialize configuration for next iteration
        iter = iter + 1;
        app.d{iter} = app.d{iter-1};
        app.p{iter} = app.p{iter-1};
        
        f=figure(6);clf(f);figHandle=axes;cla(figHandle); showHistory(app.constraints,app.classification,app.d{iter},figHandle); 
    end
end
%------------- END CODE --------------
