function [xpNames,data] = read_experiments(resultpath, varargin)

parse = inputParser;
parse.addRequired('resultpath');
parse.addOptional('runsToShow',0);
parse.parse(resultpath, varargin{:});
runsToShow = parse.Results.runsToShow;

data = {};output = {};
warning('off','all');

fsPaths = dir(resultpath);
fsPaths(~[fsPaths.isdir]) = [];
fsPaths(1:2)  = [];
xpNames = {};

if runsToShow==0; runsToShow = 1:length(fsPaths);end
for i=1:length(runsToShow)
    % fp = runsToShow(i);
    %curpath = [resultpath '/' fsPaths(fp).name];
    curpath = [resultpath '/' int2str(runsToShow(i))];    
    fsNodes = dir(curpath);
    fsNodes([fsNodes.isdir]) = []; 
    % Load results from mat files
    for f=1:length(fsNodes)
        if strcmp(fsNodes(f).name(end-2:end),'mat')
            disp(['Loading ' fsNodes(f).folder '/' fsNodes(f).name]);
            output = load([fsNodes(f).folder '/' fsNodes(f).name]);
            % Clean empty entrees for prediction and acquisition map
            nonempty_elems = arrayfun(@(s) ~isempty(s.edges),output.output.predMap);
            output.output.predMap = output.output.predMap(nonempty_elems);
            nonempty_elems = arrayfun(@(s) ~isempty(s.edges),output.output.acqMap);
            output.output.acqMap = output.output.acqMap(nonempty_elems);
            
            data{end+1} = output.output;
            name = strrep(fsNodes(f).name,'results_','');
            name = regexprep(name,'_RUN.*','');
            xpNames{end+1} = strrep(name,'_','-');
        end
    end
end
warning('on','all');

end