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


replayScoreMethods = {'rt_ts'; 'rt_ui'; 'rt_pf'; 'rt_ds'; ...
                      'wc_ts'; 'wc_ui'; 'wc_pf'; 'wc_ds'};

replayScoreMethods_fullName = {'radon Integral - wPBE time swap'; ...
                               'radon Integral - unit ID shuffle'; ...
                               'radon Integral - place field shuffle'; ...
                               'radon Integral - column cycle shuffle'; ...
                               'weighted Corr - wPBE time swap'; ...
                               'weighted Corr - unit ID shuffle'; ...
                               'weighted Corr - place field shuffle'; ...
                               'weighted Corr - column cycle shuffle'};
                           
PF_peak_min = 1;
minPBEsParticipation.pre  = 100;
minPBEsParticipation.run  = 5;
minPBEsParticipation.post = 100;



for isess = 9 %numel(sessionNames)-1

    sessionName = sessionNames{isess}
    
    sessionDir = fullfile(parentDir, sessionName);
    
%     % overall assembly tuning
%     load(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_50ms.mat']), 'assemblyTunings')
%     
%     % assembly tuning vs replay
%     load(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores_50ms.mat']), 'assemblyTunings_sub')

    
    
    % assembly tuning vs time
    
    for iEpoch = 1:3
        currEpoch = epochNames{iEpoch};
        
        s = load(fullfile(sessionDir, 'assemblyTunings', sprintf('%s.assemblyTuning_vs_time4_%s.mat', sessionName, currEpoch)));
        
        assemblyTunings_time.(currEpoch)                  = s.assemblyTunings_time;
        assemblyTuningCorrMat_time.(currEpoch)            = s.assemblyTuningCorrMat_time;
        assemblyTuningCorrMat_time_zscore.(currEpoch)     = s.assemblyTuningCorrMat_time_zscore;
        assemblyTuningPFcorr_time.(currEpoch)             = s.assemblyTuningPFcorr_time;
        assemblyTuningPFcorr_time_zscore.(currEpoch)      = s.assemblyTuningPFcorr_time_zscore;
        assemblyTunings_time_zscore.(currEpoch)           = s.assemblyTunings_time_zscore;
        assemblyTuningSpatialInfo_time.(currEpoch)        = s.assemblyTuningSpatialInfo_time;
        assemblyTuningSpatialInfo_time_zscore.(currEpoch) = s.assemblyTuningSpatialInfo_time_zscore;

        binCenters.(currEpoch) = s.binCenters;
    end

    
    
    load(fullfile(sessionDir, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
    spikes = spikes_pyr;
    epochs = fileInfo.behavior.time;
    
    if ismember(isess, 7:9)
        epochs(3,2) = epochs(3,1)+4*60*60; % the assembly tunings were calculated fo only PBEs in the first four hours
    end
    

    load(fullfile(sessionDir, 'BayesianDecoding', [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')
    PBEs = PBEInfo_replayScores;

    
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    nUnits  = numel(spikes); 

    % non-directional spatial tuning
    spatialTunings_merge = zeros(nUnits, nPosBins);
    
    for iUnit = 1:nUnits
        spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
    end    
    
    PF_prefPos = zeros(nUnits, 1);
    for iUnit = 1:nUnits
       [~, PF_prefPos(iUnit)] = max(spatialTunings_merge(iUnit, :), [], 2); 
    end

    
    
    for iepoch = 1:numel(epochNames)
        
        currEpoch = epochNames{iepoch};

        currPBEs = PBEs([PBEs.peakT] > epochs(iepoch, 1) & [PBEs.peakT] < epochs(iepoch, 2));
        
        nPBEs = numel(currPBEs);
        
        participation = zeros(nUnits, nPBEs);
        for ipbe = 1:nPBEs
            participation(:, ipbe) = sum(currPBEs(ipbe).fr_20msbin, 2) > 0;            
        end
        
        nParticipatedPBEs.(currEpoch) = sum(participation, 2);
        activeUnits.(currEpoch) = find(max(spatialTunings_merge, [], 2) > PF_peak_min & nParticipatedPBEs.(currEpoch) > minPBEsParticipation.(currEpoch));
        

        % Assembly tunings using all PBEs
        
        % real data
        [asTuning_PF_prefPosDist.(currEpoch).data, asTuning_prefPos.(currEpoch).data] = calDistPrefPosition(assemblyTunings.(currEpoch).data, PF_prefPos, activeUnits.(currEpoch));
        [asTuning_PF_prefPos_spearCorr.(currEpoch).data, asTuning_PF_prefPos_pval.(currEpoch).data] = calcorr(asTuning_prefPos.(currEpoch).data, PF_prefPos, activeUnits.(currEpoch));

        

        % surrogate data
        [asTuning_PF_prefPosDist.(currEpoch).ui, asTuning_prefPos.(currEpoch).ui] = calDistPrefPosition(assemblyTunings.(currEpoch).ui, PF_prefPos, activeUnits.(currEpoch));
        [asTuning_PF_prefPos_spearCorr.(currEpoch).ui, asTuning_PF_prefPos_pval.(currEpoch).ui] = calcorr(asTuning_prefPos.(currEpoch).ui, PF_prefPos, activeUnits.(currEpoch));


        
        % z-scores
        asTuning_PF_prefPosDist_zscore.(currEpoch) = calDistZScore(asTuning_PF_prefPosDist.(currEpoch).data, asTuning_PF_prefPosDist.(currEpoch).ui,  activeUnits.(currEpoch));
        asTuning_PF_prefPos_spearCorr_zscore.(currEpoch) = (asTuning_PF_prefPos_spearCorr.(currEpoch).data - nanmean(asTuning_PF_prefPos_spearCorr.(currEpoch).ui)) / nanstd(asTuning_PF_prefPos_spearCorr.(currEpoch).ui);
        
        
        
        
        % Assembly tunnigs calculated based on PBEs with scores in different quartiles of replay score
        
        replayScoreMethods = fieldnames(assemblyTunings_sub.pre.data);
        
        for irsm = 1:numel(replayScoreMethods)
            
            currReplayScoreMethod = replayScoreMethods{irsm};
            
            
            asTuning_PF_prefPosDist_sub.(currEpoch).(currReplayScoreMethod).data       = zeros(nUnits, 4);
            asTuning_prefPos_sub.(currEpoch).(currReplayScoreMethod).data              = zeros(nUnits, 4);
            asTuning_PF_prefPos_spearCorr_sub.(currEpoch).(currReplayScoreMethod).data = zeros(4,1);
            asTuning_PF_prefPos_pval_sub.(currEpoch).(currReplayScoreMethod).data      = zeros(4,1);
            
            asTuning_PF_prefPosDist_sub.(currEpoch).(currReplayScoreMethod).ui       = zeros(nUnits, nShuffles, 4);
            asTuning_prefPos_sub.(currEpoch).(currReplayScoreMethod).ui              = zeros(nUnits, nShuffles, 4);
            asTuning_PF_prefPos_spearCorr_sub.(currEpoch).(currReplayScoreMethod).ui = zeros(nShuffles, 4);
            asTuning_PF_prefPos_pval_sub.(currEpoch).(currReplayScoreMethod).ui      = zeros(nShuffles, 4);
            
            
            asTuning_PF_prefPosDist_sub_zscore.(currEpoch).(currReplayScoreMethod)       = zeros(nUnits, 4);
            asTuning_PF_prefPos_spearCorr_sub_zscore.(currEpoch).(currReplayScoreMethod) = zeros(4,1);
    
            
            for iPBEset = 1:4
                
                currAssemblyTuning_data    = assemblyTunings_sub.(currEpoch).data.(currReplayScoreMethod)(:, :, iPBEset);
                currAssemblyTuning_shuffle = assemblyTunings_sub.(currEpoch).ui.(currReplayScoreMethod){iPBEset};
                

                % real data
                [asTuning_PF_prefPosDist_data, asTuning_prefPos_sub.(currEpoch).(currReplayScoreMethod).data(:, iPBEset)] = calDistPrefPosition(currAssemblyTuning_data, PF_prefPos, activeUnits.(currEpoch));
                asTuning_PF_prefPosDist_sub.(currEpoch).(currReplayScoreMethod).data(:, iPBEset) = asTuning_PF_prefPosDist_data;
                
                [asTuning_PF_prefPos_spearCorr_sub.(currEpoch).(currReplayScoreMethod).data(iPBEset), asTuning_PF_prefPos_pval_sub.(currEpoch).(currReplayScoreMethod).data(iPBEset)] = calcorr(asTuning_prefPos_sub.(currEpoch).(currReplayScoreMethod).data(:, iPBEset), PF_prefPos, activeUnits.(currEpoch));

                
                % shuffle data
                [asTuning_PF_prefPosDist_ui, asTuning_prefPos_sub.(currEpoch).(currReplayScoreMethod).ui(:, :, iPBEset)] = calDistPrefPosition(currAssemblyTuning_shuffle, PF_prefPos, activeUnits.(currEpoch));
                asTuning_PF_prefPosDist_sub.(currEpoch).(currReplayScoreMethod).ui(:, :, iPBEset) = asTuning_PF_prefPosDist_ui;
                
                [asTuning_PF_prefPos_spearCorr_sub.(currEpoch).(currReplayScoreMethod).ui(:, iPBEset), asTuning_PF_prefPos_pval_sub.(currEpoch).(currReplayScoreMethod).ui(:, iPBEset)] = calcorr(asTuning_prefPos_sub.(currEpoch).(currReplayScoreMethod).ui(:, :, iPBEset), PF_prefPos, activeUnits.(currEpoch));

                
                % z-scores
                asTuning_PF_prefPosDist_sub_zscore.(currEpoch).(currReplayScoreMethod)(:, iPBEset)    = calDistZScore(asTuning_PF_prefPosDist_data, asTuning_PF_prefPosDist_ui, activeUnits.(currEpoch));
                asTuning_PF_prefPos_spearCorr_sub_zscore.(currEpoch).(currReplayScoreMethod)(iPBEset) = (asTuning_PF_prefPos_spearCorr_sub.(currEpoch).(currReplayScoreMethod).data(iPBEset) - nanmean(asTuning_PF_prefPos_spearCorr_sub.(currEpoch).(currReplayScoreMethod).ui(:, iPBEset))) / nanstd(asTuning_PF_prefPos_spearCorr_sub.(currEpoch).(currReplayScoreMethod).ui(:, iPBEset));
                
            end
        end
        

        
        % assembly tuning versus time (in 15 minutes bouts)
        
        nTimeBins = numel(binCenters.(currEpoch));
        
        asTuning_PF_prefPosDist_time.(currEpoch).data       = zeros(nUnits, nTimeBins);
        asTuning_prefPos_time.(currEpoch).data              = zeros(nUnits, nTimeBins);
        asTuning_PF_prefPos_spearCorr_time.(currEpoch).data = zeros(nTimeBins,1);
        asTuning_PF_prefPos_pval_time.(currEpoch).data      = zeros(nTimeBins,1);

        asTuning_PF_prefPosDist_time.(currEpoch).ui       = zeros(nUnits, nShuffles, nTimeBins);
        asTuning_prefPos_time.(currEpoch).ui              = zeros(nUnits, nShuffles, nTimeBins);
        asTuning_PF_prefPos_spearCorr_time.(currEpoch).ui = zeros(nShuffles, nTimeBins);
        asTuning_PF_prefPos_pval_time.(currEpoch).ui      = zeros(nShuffles, nTimeBins);
        
        asTuning_PF_prefPosDist_time_zscore.(currEpoch)       = zeros(nUnits, nTimeBins);
        asTuning_PF_prefPos_spearCorr_time_zscore.(currEpoch) = zeros(nTimeBins,1);
        
        for iBin = 1:nTimeBins
            
            currAssemblyTuning_data    = assemblyTunings_time.(currEpoch).data(:, :, iBin);
            currAssemblyTuning_shuffle = squeeze(assemblyTunings_time.(currEpoch).ui(:, :, iBin, :));
            
            
            % real data
            [asTuning_PF_prefPosDist_time.(currEpoch).data(:, iBin), asTuning_prefPos_time.(currEpoch).data(:, iBin)] = calDistPrefPosition(currAssemblyTuning_data, PF_prefPos, activeUnits.(currEpoch));
            [asTuning_PF_prefPos_spearCorr_time.(currEpoch).data(iBin), asTuning_PF_prefPos_pval_time.(currEpoch).data(iBin)] = calcorr(asTuning_prefPos_time.(currEpoch).data(:, iBin), PF_prefPos, activeUnits.(currEpoch));


            % surrogate data
            [asTuning_PF_prefPosDist_time.(currEpoch).ui(:,:,iBin), asTuning_prefPos_time.(currEpoch).ui(:, :, iBin)] = calDistPrefPosition(currAssemblyTuning_shuffle, PF_prefPos, activeUnits.(currEpoch));
            [asTuning_PF_prefPos_spearCorr_time.(currEpoch).ui(:,iBin), asTuning_PF_prefPos_pval_time.(currEpoch).ui(:,iBin)] = calcorr(asTuning_prefPos_time.(currEpoch).ui(:, :, iBin), PF_prefPos, activeUnits.(currEpoch));


            % z-scores
            asTuning_PF_prefPosDist_time_zscore.(currEpoch)(:, iBin) = calDistZScore(asTuning_PF_prefPosDist_time.(currEpoch).data(:, iBin), asTuning_PF_prefPosDist_time.(currEpoch).ui(:,:,iBin), activeUnits.(currEpoch));
            asTuning_PF_prefPos_spearCorr_time_zscore.(currEpoch)(iBin) = (asTuning_PF_prefPos_spearCorr_time.(currEpoch).data(iBin) - nanmean(asTuning_PF_prefPos_spearCorr_time.(currEpoch).ui(:,iBin))) / nanstd(asTuning_PF_prefPos_spearCorr_time.(currEpoch).ui(:,iBin));
        
        end

    end
    
%     save(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_50ms.mat']), 'asTuning_prefPos', 'asTuning_PF_prefPosDist', 'asTuning_PF_prefPosDist_zscore', 'asTuning_PF_prefPos_spearCorr', 'asTuning_PF_prefPos_pval', 'asTuning_PF_prefPos_spearCorr_zscore', 'PF_peak_min', 'minPBEsParticipation', 'activeUnits', '-append')
%     save(fullfile(sessionDir, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores_50ms.mat']), 'asTuning_prefPos_sub', 'asTuning_PF_prefPosDist_sub', 'asTuning_PF_prefPosDist_sub_zscore', 'asTuning_PF_prefPos_spearCorr_sub', 'asTuning_PF_prefPos_pval_sub', 'asTuning_PF_prefPos_spearCorr_sub_zscore', 'PF_peak_min', 'minPBEsParticipation', 'activeUnits', '-append') 
    
    save(fullfile(sessionDir, 'assemblyTunings', sprintf('%s.assemblyTuning_vs_time4.mat', sessionName)), ...
        'assemblyTunings_time', ...
        'assemblyTuningCorrMat_time', ...
        'assemblyTuningCorrMat_time_zscore', ...
        'assemblyTuningPFcorr_time', ...
        'assemblyTuningPFcorr_time_zscore', ...
        'assemblyTunings_time_zscore', ...
        'assemblyTuningSpatialInfo_time', ...
        'assemblyTuningSpatialInfo_time_zscore', ...
        'asTuning_PF_prefPosDist_time', ...
        'asTuning_prefPos_time', ...
        'asTuning_PF_prefPos_spearCorr_time', ...
        'asTuning_PF_prefPos_pval_time', ...
        'asTuning_PF_prefPosDist_time_zscore', ...
        'asTuning_PF_prefPos_spearCorr_time_zscore', ...
        'activeUnits', ...
        'binCenters', '-v7.3')
end



%% functions
        
function [prefPosDist, asTuning_prefPos] = calDistPrefPosition(assemblyTunings, PF_prefPos, activeUnits)  
    
nUnits = numel(PF_prefPos);
nSets = size(assemblyTunings, 3);
nPosBins = size(assemblyTunings, 2);

prefPosDist   = nan(nUnits, nSets);
asTuning_prefPos = nan(nUnits, nSets);

for n = 1:nSets
   currAsTunings = assemblyTunings(:,:, n);
   [~, asTuning_prefPos(:, n)] = max(currAsTunings, [], 2);
   prefPosDist(activeUnits, n) = abs(asTuning_prefPos(activeUnits, n) - PF_prefPos(activeUnits));
   prefPosDist(activeUnits, n) = prefPosDist(activeUnits, n)/nPosBins;
end

end


function asTuning_PF_prefPosDist_zscore = calDistZScore(asTuning_prefPosDist_data, asTuning_prefPosDist_ui, activeUnits)

nUnits    = size(asTuning_prefPosDist_data, 1);
% nShuffles = size(asTuning_prefPosDist_ui, 2);

nActiveUnits = numel(activeUnits);
% asTuning_PF_prefPosDist_prctl = nan(nUnits, 1);
asTuning_PF_prefPosDist_zscore = nan(nUnits, 1);
for iUnit = 1:nActiveUnits
    
    currUnit = activeUnits(iUnit);
%     asTuning_PF_prefPosDist_prctl(currUnit) = numel(find(asTuning_prefPosDist_data(currUnit) < asTuning_prefPosDist_ui(currUnit, :)))/ nShuffles;
    asTuning_PF_prefPosDist_zscore(currUnit) = (asTuning_prefPosDist_data(currUnit)- nanmean(asTuning_prefPosDist_ui(currUnit, :)))/ nanstd(asTuning_prefPosDist_ui(currUnit, :));

end

end
  
function [corrCoeff, pval] = calcorr(asTuning_prefPos, PF_prefPos, activeUnits)

nSets = size(asTuning_prefPos, 2);

corrCoeff = nan(nSets, 1);
pval      = nan(nSets, 1);

for n = 1:nSets   
    [corrCoeff(n), pval(n)] = corr(asTuning_prefPos(activeUnits, n), PF_prefPos(activeUnits), 'type', 'spearman');
end

end
    
