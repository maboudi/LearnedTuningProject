function assemblyTuning_continuous(sessionNumber)


% we might need to consider only stable units when calculating the assembly
% tuning for a specific unit

% sz = getenv('SLURM_CPUS_PER_TASK');
% p  = parpool('local', str2double(sz));


% addpath(genpath('/nfs/turbo/umms-kdiba/NCMLproject/ReplayPreplayAnalyses'))


% parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/concat_GreatLakes_datasets_temp/';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name;


basePath    = fullfile(parentDir, sessionName);
storagePath = fullfile(basePath, 'assemblyTunings');
mkdir(storagePath)


load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')


% behavioral epochs
behavior = fileInfo.behavior;

startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
startT.post = behavior.time(2,2); endT.post = behavior.time(2,2)+3*60*60; 


clusterQuality = clusterQuality; 


spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

nUnits  = numel(spikes); 


% non-directional spatial tuning
spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
end


%% assembly Tunings considering all of the PBEs
% Spatial information and correlation with place fields of assesmbly tunings
% and corresponding z-scores by comparing the values againts distibutions calculated based on unit ID shuffle surrogates


periodNames = {'pre'; 'run'; 'post'};

winDur    = 900;
stepDur   = 300;

nUIshuffles = 10;
nSeg        = 1;

nUIshuffles_perSeg = nUIshuffles/nSeg;

for iperiod = 3
    
    currPeriod = periodNames{iperiod};
    epoch     = [startT.(currPeriod) endT.(currPeriod)];

    fprintf(['/nProcessing ' currPeriod ' ..'])
    
    
    
    % actual data
    fprintf('\nActual PBEs')
    
    ifshuffle = 0;
    [currAssemblyTunings, px, winCenters.(currPeriod)] = calculateAssemblyTuning_continuous(spikes, epoch, winDur, stepDur, clusterQuality, ifshuffle);
    
    virtualOccupancy.(currPeriod) = px/sum(px); %%% correct this; px should be recalculated using all units
    nWins = numel(winCenters.(currPeriod));
    
    assemblyTunings_contin.(currPeriod)  = currAssemblyTunings;
    
    assemblyTuningCorrMat_contin.(currPeriod)     = zeros(nUnits, nUnits, nWins);
    assemblyTuningPFcorr_contin.(currPeriod)      = zeros(nUnits, nWins);
    assemblyTuningSpatialInfo_contin.(currPeriod) = zeros(nUnits, nWins);   

    for iwin = 1:nWins
        
        assemblyTuningCorrMat_contin.(currPeriod)(:, :, iwin) = corr(currAssemblyTunings(:, :, iwin)', spatialTunings_merge'); 
        assemblyTuningPFcorr_contin.(currPeriod)(:, iwin)     = diag(assemblyTuningCorrMat_contin.(currPeriod)(:, :, iwin)); % the correlation of PF and assembly tuning for a given unit

        meanFR = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* currAssemblyTunings(:, :, iwin), 2) ./ sum(virtualOccupancy.(currPeriod)); % The calculation of mean firing rate could be a bit questionable
        assemblyTuningSpatialInfo_contin.(currPeriod)(:, iwin) = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* (currAssemblyTunings(:,:,iwin) ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings(:,:,iwin) ./ repmat(meanFR, [1 nPosBins])), 2);
    end
    
    
    % surrogate (unit Id shuffle)
     fprintf('\nUnit ID shuffle PBEs')
    

    ifshuffle = 1;
    
    currAssemblyTunings = zeros(nUnits, nPosBins, nWins, nUIshuffles); % we are not going to store this variable, only need it for z-scoring 
    currCorrMat         = zeros(nUnits, nUnits, nWins, nUIshuffles);
    currPFcorr          = zeros(nUnits, nWins, nUIshuffles);
    currSpatialInfo     = zeros(nUnits, nWins, nUIshuffles);
    currOccupancy       = virtualOccupancy.(currPeriod);

    
    for seg = 1:nSeg
        fprintf('\nSegment %d out of %d', seg, nSeg)
        parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)

            shuffleAssemblyTunings = calculateAssemblyTuning_continuous(spikes, epoch, winDur, stepDur, clusterQuality, ifshuffle);
            currAssemblyTunings(:, :, :, i_ui)= shuffleAssemblyTunings;
            
            for iwin = 1:nWins
                currCorrMat(:,:, iwin, i_ui) = corr(shuffleAssemblyTunings(:,:,iwin)', spatialTunings_merge'); 
                currPFcorr(:, iwin, i_ui)    = diag(currCorrMat(:,:, iwin, i_ui));

                meanFR                         = sum(repmat(currOccupancy', [nUnits 1]).* shuffleAssemblyTunings(:, :, iwin), 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
                currSpatialInfo(:, iwin, i_ui) = sum(repmat(currOccupancy', [nUnits 1]).* (shuffleAssemblyTunings(:,:,iwin) ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings(:,:,iwin) ./ repmat(meanFR, [1 nPosBins])), 2);
            end
        end
    end

    
    % comparing actual assembly tunings with the assembly tunings based on
    % unit ID shuffles
    assemblyTunings_contin_zscore.(currPeriod)  = (assemblyTunings_contin.(currPeriod) - mean(currAssemblyTunings, 4)) ./ std(currAssemblyTunings, [], 4);
    assemblyTuningCorrMat_contin_zscore.(currPeriod)   = (assemblyTuningCorrMat_contin.(currPeriod) - mean(currCorrMat, 4)) ./ std(currCorrMat, [], 4);
    
    assemblyTuningPFcorr_contin_zscore.(currPeriod)      = (assemblyTuningPFcorr_contin.(currPeriod) - mean(currPFcorr, 3)) ./ std(currPFcorr, [], 3);
    assemblyTuningSpatialInfo_contin_zscore.(currPeriod) = (assemblyTuningSpatialInfo_contin.(currPeriod) - mean(currSpatialInfo, 3)) ./ std(currSpatialInfo, [], 3);
    
    
    fileName = [sessionName '.assemblyTuning_vs_contin.mat'];
    save(fullfile(storagePath, fileName),...
    'assemblyTunings_contin', 'assemblyTuningCorrMat_contin', 'assemblyTuningPFcorr_contin', 'assemblyTuningSpatialInfo_contin', ...
    'assemblyTunings_contin_zscore', 'assemblyTuningCorrMat_contin_zscore', 'assemblyTuningPFcorr_contin_zscore', 'assemblyTuningSpatialInfo_contin_zscore', 'winCenters')

    
end



end

