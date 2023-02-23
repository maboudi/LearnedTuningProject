function BayesianDecoding_June21(sessionNumber) 

addpath(genpath('/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/ReplayPreplayAnalyses'))

sz = getenv('SLURM_CPUS_PER_TASK');
p = parpool('local', str2num(sz));


% parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/sessions_calculated_PBEinfo3';

parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/sessions_calculated_PBEinfo3';


files = dir(parentDir);
isdir = [files.isdir];
subfolders = files(isdir);

sessionNames = subfolders(3:end);



sessionName = sessionNames(sessionNumber).name;


%%  Bayesian decoding

basePath = fullfile(parentDir, sessionName);


% spike info
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']))
if exist('spikes_pyr', 'var')
    clear spikes
    spikes = spikes_pyr;
end




% PBEs
load(fullfile(basePath, 'PopulationBurstEvents', [sessionName '.PBEInfo.mat']))



binDur = 0.02;
[PBEInfo_replayScores, PBEInfo_null, sigMatrix] = BDreplayDetect_mixedDir(PBEInfo_Bayesian, spikes, fileInfo, binDur);



storagePath = fullfile(basePath, 'BayesianDecoding');
if ~exist(storagePath, 'dir'); mkdir(storagePath); end

save(fullfile(storagePath, [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores', 'PBEInfo_null', 'sigMatrix', '-v7.3')

