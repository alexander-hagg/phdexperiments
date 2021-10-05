setSchedulerMessageHandler(@disp);

pc = parcluster('local');
% explicitly set the JobStorageLocation to the temp directory that was created in your sbatch script
mkdir('/tmp/testMatlabtmp');
pc.JobStorageLocation = '/tmp/testMatlabtmp';
parpool(pc, 30);

