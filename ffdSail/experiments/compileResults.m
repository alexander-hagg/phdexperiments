clear
rep = {'Parsec','FFD','CPPN'};
for iRep=1:length(rep)
%% Unpack and get medians
%
%   predicted:
%       fitness
%       drag
%       lift        
%       
%   true:
%       fitness
%       drag
%       lift
%
%   median:
%       parameter values
%       performing designs
%       compared to CMA-ES Parsec
%

%% Compile Values
load(['~/Code/data/ffdSail/sail' rep{iRep} '_true.mat']);
allOutput = output;
clearvars predFitness predcD predcL trueFitness truecD truecL genes
for iRun = 1:length(allOutput)
    output = allOutput{iRun}; eval(output.unpack)
    evalPt = p.nInitialSamples:p.data.mapEvalMod:p.nTotalSamples;
    for iItr = 1:length(evalPt)
        predFitness(:,:,iItr,iRun)  = predMap(evalPt(iItr)).fitness;
        predcD(:,:,iItr,iRun)       = predMap(evalPt(iItr)).cD;
        predcL(:,:,iItr,iRun)       = predMap(evalPt(iItr)).cL;
        
        trueFitness(:,:,iItr,iRun)  = predMap(evalPt(iItr)).fitness_true;
        truecD(:,:,iItr,iRun)       = predMap(evalPt(iItr)).cD_true;
        truecL(:,:,iItr,iRun)       = predMap(evalPt(iItr)).cL_true;
        
        genes(:,:,:,iItr,iRun)      = predMap(evalPt(iItr)).genes;
    end
end

% Predicted Values
eval(['all' rep{iRep} '.predFitness = predFitness;']);
eval(['all' rep{iRep} '.predcD      = predcD;']);
eval(['all' rep{iRep} '.predcL      = predcL;']);

% True Values
eval(['all' rep{iRep} '.trueFitness = trueFitness;']);
eval(['all' rep{iRep} '.truecD      = truecD;']);
eval(['all' rep{iRep} '.truecL      = truecL;']);

% Errors (absolute)
eval(['all' rep{iRep} '.absErrorFit = ' 'all' rep{iRep} '.predFitness - ' 'all' rep{iRep} '.trueFitness'])
eval(['all' rep{iRep} '.absErrorcD  = ' 'all' rep{iRep} '.predcD      - ' 'all' rep{iRep} '.truecD'     ])
eval(['all' rep{iRep} '.absErrorcL  = ' 'all' rep{iRep} '.predcL      - ' 'all' rep{iRep} '.truecL'     ])

% Errors (%)
eval(['all' rep{iRep} '.percErrorFit = abs(100*(abs(' 'all' rep{iRep} '.predFitness - ' 'all' rep{iRep} '.trueFitness)./' 'all' rep{iRep} '.trueFitness))'])
eval(['all' rep{iRep} '.percErrorcD  = abs(100*(abs(' 'all' rep{iRep} '.predcD      - ' 'all' rep{iRep} '.truecD)     ./' 'all' rep{iRep} '.truecD))'])
eval(['all' rep{iRep} '.percErrorcL  = abs(100*(abs(' 'all' rep{iRep} '.predcL      - ' 'all' rep{iRep} '.truecL)     ./' 'all' rep{iRep} '.truecL))'])

eval(['all' rep{iRep} '.genes   = genes;']);

%% Compute Medians
% Median Values
fields = fieldnames(allParsec);
for iField = 1:length(fields)-2
    eval([rep{iRep} '.' fields{iField} ' = nanmedian(all' rep{iRep} '.' fields{iField} ',4)'])
end
eval([rep{iRep} '.trueFitStd  = nanstd(all' rep{iRep} '.trueFitness,[],4)'])

% Median and variance per iteration compared to 'optimal' parsec by cma-es
load('c3Result.mat'); optFit = cmaOpt.fitness;
for iItr = 1:length(evalPt)
      eval(['byBinPerf = ' rep{iRep} '.trueFitness(:,:,iItr)./optFit'])
      eval(['byBinStd  = ' rep{iRep}  '.trueFitStd(:,:,iItr)./optFit']);
      eval([rep{iRep} '.trueFitOptMedian(iItr) = nanmedian(byBinPerf(:))']);     
      eval([rep{iRep} '.trueFitOptStd(iItr)    = nanmedian(byBinStd (:))']);     
end

eval([rep{iRep} ' = orderfields(' rep{iRep} ');']);

% Find median performing individuals in each bin (super slow for now...)
for i=1:d.featureRes(1)
    for j=1:d.featureRes(2)
        for iItr = 1:length(evalPt)
            eval([ 'A = ' rep{iRep} '.trueFitness(i,j,iItr);'])
            eval([ 'AA = squeeze(all' rep{iRep} '.trueFitness(i,j,iItr,:));'])
            [val,idx]=min(abs((A-AA)));
            eval([ 'medianGenes = squeeze(all' rep{iRep} '.genes(i,j,:,iItr,idx));']);
            eval([rep{iRep} '.genes(i,j,:,iItr) = medianGenes;']);
        end
    end
end

eval([rep{iRep} '.evalPt  = evalPt;']);
    
save([rep{iRep} '_result.mat'], rep{iRep})
end

% plot([Parsec.trueFitOptMedian; Parsec.trueFitOptMedian-Parsec.trueFitOptStd;Parsec.trueFitOptMedian+Parsec.trueFitOptStd]','b-'); hold on;
% plot([FFD.trueFitOptMedian; FFD.trueFitOptMedian-FFD.trueFitOptStd;FFD.trueFitOptMedian+FFD.trueFitOptStd]','g-')
% plot([CPPN.trueFitOptMedian; CPPN.trueFitOptMedian-CPPN.trueFitOptStd;CPPN.trueFitOptMedian+CPPN.trueFitOptStd]','r-')
%     


