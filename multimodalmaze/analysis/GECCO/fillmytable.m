load('/scratch/ahagg2s/GECCO2019/planner_MUTATION_4096/analysis.mat');

mean(mCor(:,6:end)')
std(mCor(:,6:end)')

mean(mInCor(:,6:end)')
std(mInCor(:,6:end)')

mean(mDrift(:,6:end)')
std(mDrift(:,6:end)')

%%
load('/scratch/ahagg2s/GECCO2019/controller_MUTATION_4096/analysis.mat');
sel = 1:4;

mean(mCor(:,sel)')
std(mCor(:,sel)')

mean(mInCor(:,sel)')
std(mInCor(:,sel)')

mean(mDrift(:,sel)')
std(mDrift(:,sel)')

