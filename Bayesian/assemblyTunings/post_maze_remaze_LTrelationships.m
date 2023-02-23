
% to calcualte the relationship between learned tunings calculated during
% run, post, and re-maze


clear;
clc;


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject';

rr = dir(fullfile(parentDir, 'assemblyTuning_finalResults'));



for sessionNumber = 8

sessionName = rr(sessionNumber+2).name

basePath = fullfile(parentDir, 'assemblyTuning_finalResults', sessionName);


gw = gausswindow(3,9);

% place fields
s = load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'fileInfo', 'spikes_pyr');
fileInfo = s.fileInfo;
spikes   = s.spikes_pyr;

nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
nUnits   = numel(spikes); 


% non-directional spatial tuning

spatialTunings = nan(nUnits, nPosBins);
peakPFfiring   = nan(nUnits, 1);
peakLoc        = nan(nUnits, 1);

for iUnit = 1:nUnits
    spatialTunings(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni; 
    [~, peakLoc(iUnit)] = max(spatialTunings(iUnit, :));
    
    peakPFfiring(iUnit) = spikes(iUnit).peakFR.uni;

end
[~, PFSortIdx] = sort(peakLoc, 'ascend');



% % % learned tunings

% maze
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_maze.mat']), 'learnedTunings', 'learnedTuningPFcorr');

maze_LTs      = s.learnedTunings;
maze_LTPFcorr = s.learnedTuningPFcorr;


% for iUnit = 1:nUnits
%     maze_LTs(iUnit, :) = conv(maze_LTs(iUnit, :), gw, 'same');
% end


clear s
% post

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), 'assemblyTunings', 'assemblyTuningPFcorr');

post_LTs      = s.assemblyTunings.post.data;
post_LTPFcorr = s.assemblyTuningPFcorr.post.data;

for iUnit = 1:nUnits
    post_LTs(iUnit, :) = conv(post_LTs(iUnit, :), gw, 'same');
end



clear s
% remaze
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_remaze.mat']), 'learnedTunings', 'learnedTuningPFcorr');

remaze_LTs      = s.learnedTunings;
remaze_LTPFcorr = s.learnedTuningPFcorr;

clear s


% corrrelation bw/ learned tuning in post and learned tuning in Maze or
% Re-maze

allCorrs       = corr(post_LTs', maze_LTs');
post_maze_corr = diag(allCorrs);


allCorrs       = corr(post_LTs', remaze_LTs');
post_remaze_corr = diag(allCorrs);





end