clear
clc

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults';

files = dir(parentDir);
isdir = [files.isdir];
subfolders = files(isdir);

sessionNames = {subfolders(3:end).name};

epochNames = {'pre'; 'run'; 'post'};
posBinSize = 2; % in cm
nShuffles  = 100; % number of surrogate assembly tunings calculated based on Unit-ID shuffle


replayScoreMethods = {'rt_ui'; 'wc_ui'; 'wc_ts'};


replayScoreMethods_fullName = {'radon Integral - unit ID shuffle'; ...
                               'weighted Corr - unit ID shuffle' ; ...
                               'weighted Corr - wPBE time swap'};
                           


for isess = [6:15]%1:numel(sessionNames)-1

    sessionName = sessionNames{isess}
    
    sessionDir = fullfile(parentDir, sessionName);
    
    % overall assembly tuning
    load(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'assemblyTunings')
    
    % assembly tuning vs replay
    load(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores.mat']), 'assemblyTunings_sub')

    
    
    % assembly tuning vs time
    
    load(fullfile(sessionDir, 'assemblyTunings', sprintf('%s.assemblyTuning_vs_time4.mat', sessionName)), 'assemblyTunings_time', 'binCenters');

    
    
    load(fullfile(sessionDir, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
    spikes = spikes_pyr;
    
    
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    nUnits  = numel(spikes); 

    % non-directional spatial tuning
    spatialTunings = zeros(nUnits, nPosBins);
    PF_prefPos = zeros(nUnits, 1);

    for iUnit = 1:nUnits
        spatialTunings(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
        [~, PF_prefPos(iUnit)] = max(spatialTunings(iUnit, :), [], 2); 
    end    
    
    
    for iepoch = 1:numel(epochNames)
        
        currEpoch = epochNames{iepoch};


        % Assembly tunings using all PBEs
        
        
        asTuning_PF_KLdiv.(currEpoch).data   = zeros(nUnits, 1);
        asTuning_PF_KLdiv.(currEpoch).ui     = zeros(nUnits, nShuffles);
        asTuning_PF_KLdiv_zscore.(currEpoch) = zeros(nUnits, 1);
        for iUnit = 1:nUnits
            
            % real data
            asTuning_PF_KLdiv.(currEpoch).data(iUnit) = calKLDivergence(assemblyTunings.(currEpoch).data(iUnit, :), spatialTunings(iUnit, :));
            
            
            % surrogate data
            currLearnedTunings = assemblyTunings.(currEpoch).ui(iUnit, :, :);
            currLearnedTunings = permute(currLearnedTunings, [2 3 1]);
            
            asTuning_PF_KLdiv.(currEpoch).ui(iUnit, :) = calKLDivergence(currLearnedTunings, spatialTunings(iUnit, :));
            
        end
        
        asTuning_PF_KLdiv_zscore.(currEpoch) = (asTuning_PF_KLdiv.(currEpoch).data - mean(asTuning_PF_KLdiv.(currEpoch).ui, 2))./std(asTuning_PF_KLdiv.(currEpoch).ui ,[], 2);
            
        
        
        
        
        % Assembly tunnigs calculated based on PBEs with scores in different quartiles of replay score
        
        replayScoreMethods = fieldnames(assemblyTunings_sub.pre.data);
        
        for irsm = 1:numel(replayScoreMethods)
            
            currReplayScoreMethod = replayScoreMethods{irsm};
            
            asTuning_PF_KLdiv_sub.(currEpoch).(currReplayScoreMethod).data   = zeros(nUnits, 4);
            asTuning_PF_KLdiv_sub.(currEpoch).(currReplayScoreMethod).ui     = zeros(nUnits, nShuffles, 4);
            
            asTuning_PF_KLdiv_sub_zscore.(currEpoch).(currReplayScoreMethod) = zeros(nUnits, 4);
            
            
            for iPBEset = 1:4
                
                currAssemblyTuning_data    = assemblyTunings_sub.(currEpoch).data.(currReplayScoreMethod)(:, :, iPBEset);
                currAssemblyTuning_shuffle = assemblyTunings_sub.(currEpoch).ui.(currReplayScoreMethod){iPBEset};
                
                
                for iUnit = 1:nUnits
                    
                    % real data
                    asTuning_PF_KLdiv_sub.(currEpoch).(currReplayScoreMethod).data(iUnit, iPBEset) = calKLDivergence(currAssemblyTuning_data(iUnit, :), spatialTunings(iUnit, :));
                    
                    % surrogate data
                    currLearnedTunings = currAssemblyTuning_shuffle(iUnit, :, :);
                    currLearnedTunings = permute(currLearnedTunings, [2 3 1]);
                    
                    asTuning_PF_KLdiv_sub.(currEpoch).(currReplayScoreMethod).ui(iUnit, :, iPBEset) = calKLDivergence(currLearnedTunings, spatialTunings(iUnit, :));
                    
                end
                
                asTuning_PF_KLdiv_sub_zscore.(currEpoch).(currReplayScoreMethod)(:, iPBEset) = (asTuning_PF_KLdiv_sub.(currEpoch).(currReplayScoreMethod).data(:, iPBEset) ...
                    - mean(asTuning_PF_KLdiv_sub.(currEpoch).(currReplayScoreMethod).ui(:, :, iPBEset), 2))./ std(asTuning_PF_KLdiv_sub.(currEpoch).(currReplayScoreMethod).ui(:, :, iPBEset), [], 2);

            end
        end
        

        
        % assembly tuning versus time (in 15 minutes bouts)
        
        nTimeBins = numel(binCenters.(currEpoch));
        
        asTuning_PF_KLdiv_time.(currEpoch).data   = zeros(nUnits, nTimeBins);
        asTuning_PF_KLdiv_time.(currEpoch).ui     = zeros(nUnits, nShuffles, nTimeBins);
        asTuning_PF_KLdiv_time_zscore.(currEpoch) = zeros(nUnits, nTimeBins);
        
        
        for iBin = 1:nTimeBins
            
            currAssemblyTuning_data    = assemblyTunings_time.(currEpoch).data(:, :, iBin);
            currAssemblyTuning_shuffle = squeeze(assemblyTunings_time.(currEpoch).ui(:, :, iBin, :));
            
            
            for iUnit = 1:nUnits
                    
                % real data
                asTuning_PF_KLdiv_time.(currEpoch).data(iUnit, iBin) = calKLDivergence(currAssemblyTuning_data(iUnit, :), spatialTunings(iUnit, :));

                % surrogate data
                currLearnedTunings = currAssemblyTuning_shuffle(iUnit, :, :);
                currLearnedTunings = permute(currLearnedTunings, [2 3 1]);

                asTuning_PF_KLdiv_time.(currEpoch).ui(iUnit, :, iBin) = calKLDivergence(currLearnedTunings, spatialTunings(iUnit, :));

            end

            asTuning_PF_KLdiv_time_zscore.(currEpoch)(:, iBin) = (asTuning_PF_KLdiv_time.(currEpoch).data(:, iBin) ...
                - mean(asTuning_PF_KLdiv_time.(currEpoch).ui(:, :, iBin), 2))./ std(asTuning_PF_KLdiv_time.(currEpoch).ui(:, :, iBin), [], 2);
        end

    end
    
    save(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'asTuning_PF_KLdiv', 'asTuning_PF_KLdiv_zscore', '-append')
    save(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores.mat']), 'asTuning_PF_KLdiv_sub', 'asTuning_PF_KLdiv_sub_zscore', '-append') 
    save(fullfile(sessionDir, 'assemblyTunings', sprintf('%s.assemblyTuning_vs_time4.mat', sessionName)), 'asTuning_PF_KLdiv_time', 'asTuning_PF_KLdiv_time_zscore', '-append')
    
end



%% functions

function KLdiv = calKLDivergence(learnedTunings, spatialTunings)

% this function works well when when the learnedTunings matrix contains multiple instances of learned tunings (for examples if they belong to different time windows or different instances of unit identity shuffle)

if size(spatialTunings, 2) > size(spatialTunings, 1)
    spatialTunings = spatialTunings';
end
spatialTunings = spatialTunings ./ repmat(sum(spatialTunings, 1), [size(spatialTunings, 1) 1]);
spatialTunings = spatialTunings + eps;


nPosBins = size(spatialTunings, 1);

if size(learnedTunings, 1) ~= nPosBins
    learnedTunings = learnedTunings';
end

learnedTunings = learnedTunings ./ repmat(sum(learnedTunings, 1), [size(learnedTunings, 1) 1]);
learnedTunings = learnedTunings + eps;


nLTs = size(learnedTunings, 2); % number of learned tunings

spatialTuningTerm = repmat(spatialTunings, [1 nLTs]);
KLdiv = sum(spatialTuningTerm .* (log(spatialTuningTerm) - log(learnedTunings)), 1);

end
