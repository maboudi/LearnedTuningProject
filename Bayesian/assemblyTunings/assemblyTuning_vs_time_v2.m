function assemblyTuning_vs_time_v2(sessionNumber)


% we might need to consider only stable units when calculating the assembly
% tuning for a specific unit

sz = getenv('SLURM_CPUS_PER_TASK');
p  = parpool(str2double(sz));


addpath(genpath('/nfs/turbo/umms-kdiba/NCMLproject/ReplayPreplayAnalyses'))


parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
% parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/concat_GreatLakes_datasets_temp/';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name;


basePath    = fullfile(parentDir, sessionName);
storagePath = fullfile(basePath, 'assemblyTunings');
mkdir(storagePath)


load(fullfile(basePath, 'BayesianDecoding', [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')


% behavioral epochs
behavior = fileInfo.behavior;

startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
startT.post = behavior.time(2,2); endT.post = behavior.time(3,2); 


clusterQuality = clusterQuality; 


spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);


nUnits  = numel(spikes); 



% non-directional spatial tuning
spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
end


% we are going to analyze assembly tuning separately for each epoch, so we
% need to reorganize the PBEs


% PBEInfo = PBEInfo(acceptedIdx); % only PBEs qulaified in terms of duration, number of participant units, the brain satte during which they occurred
PBEInfo = PBEInfo_replayScores;
PBEInfo = PBEInfo(ismember({PBEInfo.brainState}, {'QW'; 'NREM'}));

peakT   = [PBEInfo.peakT]';


%% assembly Tunings considering all of the PBEs
% Spatial information and correlation with place fields of assesmbly tunings
% and corresponding z-scores by comparing the values againts distibutions calculated based on unit ID shuffle surrogates


periodNames = {'pre'; 'run'; 'post'};

binDur  = 900; % 15 minutes
stepDur = 300; % 5 minutes


nUIshuffles = 100;
nSeg        = 4;

nUIshuffles_perSeg = nUIshuffles/nSeg;

for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};

    fprintf(['/nProcessing ' currPeriod ' ..'])
    
    
    allPBEs = PBEInfo(strcmp({PBEInfo.epoch}, currPeriod));
    TN_PBEs = numel(allPBEs);
    
    posteriorProbMatrix = cell(TN_PBEs, 1);
    for pbe = 1:TN_PBEs
        posteriorProbMatrix{pbe} = allPBEs(pbe).posteriorProbMat;
    end
    posteriorProbMatrix = cell2mat(posteriorProbMatrix');

    virtualOccupancy = mean(posteriorProbMatrix, 2);
    virtualOccupancy = virtualOccupancy / sum(virtualOccupancy); % normalize


    
    totalT = endT.(currPeriod) - startT.(currPeriod);
    nBins  = floor((totalT - binDur)/stepDur) + 1;

    binStarts   = startT.(currPeriod) + (0:nBins-1)*stepDur;
    binEnds     = binStarts + binDur;
    binCenters = binStarts + binDur/2;
    
    
    selectPBEs = cell(nBins, 1);
    nPBEs = zeros(nBins, 1);
    for ibin = 1:nBins
        selectPBEs{ibin} = find(peakT >= binStarts(ibin) & peakT < binEnds(ibin));
        nPBEs(ibin) = numel(selectPBEs{ibin});
    end
    
    
    
    % actual PBEs
    fprintf('\nActual PBEs')
    
    assemblyTuningCorrMat_time     = zeros(nUnits, nUnits, nBins);
    assemblyTuningPFcorr_time      = zeros(nUnits, nBins);
    assemblyTuningSpatialInfo_time = zeros(nUnits, nBins);   

    ifshuffle = 0;
    currAssemblyTunings = calculateAssemblyTuning_vs_time(PBEInfo, spikes, selectPBEs, clusterQuality, ifshuffle);
    
    
    assemblyTunings_time.data  = currAssemblyTunings;
    for ibin = 1:nBins
        
        assemblyTuningCorrMat_time(:, :, ibin) = corr(currAssemblyTunings(:, :, ibin)', spatialTunings_merge'); 
        assemblyTuningPFcorr_time(:, ibin)     = diag(assemblyTuningCorrMat_time(:, :, ibin)); % the correlation of PF and assembly tuning for a given unit

        meanFR = sum(repmat(virtualOccupancy', [nUnits 1]).* currAssemblyTunings(:, :, ibin), 2) ./ sum(virtualOccupancy); % The calculation of mean firing rate could be a bit questionable
        assemblyTuningSpatialInfo_time(:, ibin) = sum(repmat(virtualOccupancy', [nUnits 1]).* (currAssemblyTunings(:,:,ibin) ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings(:,:,ibin) ./ repmat(meanFR, [1 nPosBins])), 2);
    end
    
    
    
    % surrogate (unit Id shuffle)
    fprintf('\nUnit ID shuffle PBEs')
    

    ifShuffle = 1;
    currAssemblyTunings = zeros(nUnits, nPosBins, nBins, nUIshuffles); % we are not going to store this variable, only need it for z-scoring 
    currCorrMat         = zeros(nUnits, nUnits, nBins, nUIshuffles);
    currPFcorr          = zeros(nUnits, nBins, nUIshuffles);
    currSpatialInfo     = zeros(nUnits, nBins, nUIshuffles);
    currOccupancy       = virtualOccupancy;

    
    for seg = 1:nSeg
        fprintf('\nSegment %d out of %d', seg, nSeg)
        parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)

            shuffleAssemblyTunings = calculateAssemblyTuning_vs_time(PBEInfo, spikes, selectPBEs, clusterQuality, ifShuffle);
            currAssemblyTunings(:, :, :, i_ui)= shuffleAssemblyTunings;
            
            for ibin = 1:nBins
                currCorrMat(:,:, ibin, i_ui) = corr(shuffleAssemblyTunings(:,:,ibin)', spatialTunings_merge'); 
                currPFcorr(:, ibin, i_ui)    = diag(currCorrMat(:,:, ibin, i_ui));

                meanFR                         = sum(repmat(currOccupancy', [nUnits 1]).* shuffleAssemblyTunings(:, :, ibin), 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
                currSpatialInfo(:, ibin, i_ui) = sum(repmat(currOccupancy', [nUnits 1]).* (shuffleAssemblyTunings(:,:,ibin) ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings(:,:,ibin) ./ repmat(meanFR, [1 nPosBins])), 2);
            end
        end
    end
    
    assemblyTunings_time.ui = currAssemblyTunings;
    
    % comparing actual assembly tunings with the assembly tunings based on
    % unit ID shuffles
    assemblyTunings_time_zscore         = (assemblyTunings_time.data - mean(assemblyTunings_time.ui, 4)) ./ std(assemblyTunings_time.ui, [], 4);
    assemblyTuningCorrMat_time_zscore   = (assemblyTuningCorrMat_time - mean(currCorrMat, 4)) ./ std(currCorrMat, [], 4);
    
    assemblyTuningPFcorr_time_zscore      = (assemblyTuningPFcorr_time - mean(currPFcorr, 3)) ./ std(currPFcorr, [], 3);
    assemblyTuningSpatialInfo_time_zscore = (assemblyTuningSpatialInfo_time - mean(currSpatialInfo, 3)) ./ std(currSpatialInfo, [], 3);
    
    
    fileName = [sessionName sprintf('.assemblyTuning_vs_time4_%s.mat', currPeriod)];
    
    save(fullfile(storagePath, fileName), 'nPBEs', 'binCenters', ...
            'assemblyTunings_time', 'assemblyTuningCorrMat_time', 'assemblyTuningPFcorr_time', 'assemblyTuningSpatialInfo_time', ...
            'assemblyTunings_time_zscore', 'assemblyTuningCorrMat_time_zscore', 'assemblyTuningPFcorr_time_zscore', 'assemblyTuningSpatialInfo_time_zscore', '-v7.3')

end


end

