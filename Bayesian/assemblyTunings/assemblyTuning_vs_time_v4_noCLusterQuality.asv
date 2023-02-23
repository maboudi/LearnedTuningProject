function assemblyTuning_vs_time_v4_noCLusterQuality(sessionNumber)


% sz = getenv('SLURM_CPUS_PER_TASK');
% 
% theCluster = parcluster('local');
% 
% JobFolder = sprintf('/home/kmaboudi/.matlab/trash/job_%s', sessionNumber);
% mkdir(JobFolder)
% theCluster.JobStorageLocation = JobFolder;
% 
% p = parpool(theCluster, str2double(sz));



% addpath(genpath('/nfs/turbo/umms-kdiba/NCMLproject/ReplayPreplayAnalyses'))

% parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name


basePath    = fullfile(parentDir, sessionName);
storagePath = fullfile(basePath, 'assemblyTunings');
mkdir(storagePath)




load(fullfile(basePath, 'PopulationBurstEvents', [sessionName '.PBEInfo.mat']), 'PBEInfo_Bayesian')
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
% load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')


% behavioral epochs
behavior = fileInfo.behavior;

startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
startT.post = behavior.time(2,2); endT.post = behavior.time(3,2); 


% clusterQuality = clusterQuality; 
clusterQuality = [];

spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);


nUnits  = numel(spikes); 


% non-directional spatial tuning
spatialTunings_merge = nan(nUnits, nPosBins);
PF_prefPos = nan(nUnits, 1);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
    [~, PF_prefPos(iUnit)] = max(spatialTunings_merge(iUnit, :), [], 2);
end


% we are going to analyze assembly tuning separately for each epoch, so we
% need to reorganize the PBEs


% PBEInfo = PBEInfo(acceptedIdx); % only PBEs qulaified in terms of duration, number of participant units, the brain satte during which they occurred
PBEInfo = PBEInfo_Bayesian;
PBEInfo = PBEInfo(ismember({PBEInfo.brainState}, {'QW'; 'NREM'}));

peakT   = [PBEInfo.peakT]';


%% assembly Tunings considering all of the PBEs
% Spatial information and correlation with place fields of assesmbly tunings
% and corresponding z-scores by comparing the values againts distibutions calculated based on unit ID shuffle surrogates


epochNames = {'pre'; 'run'; 'post'};

binDur  = 900; % 15 minutes
stepDur = 300; % 5 minutes


for iepoch = 1:3
    
    currEpoch = epochNames{iepoch};
    fprintf(['\nProcessing ' currEpoch ' ..'])
    
    
    allPBEs = PBEInfo(strcmp({PBEInfo.epoch}, currEpoch));
    tNum_PBEs = numel(allPBEs);
     
%     posteriorProbMatrix = cell(tNum_PBEs, 1);
%     for pbe = 1:tNum_PBEs
%         posteriorProbMatrix{pbe} = allPBEs(pbe).posteriorProbMat;
%     end
%     posteriorProbMatrix = cell2mat(posteriorProbMatrix');
% 
%     virtualOccupancy = mean(posteriorProbMatrix, 2);
%     virtualOccupancy = virtualOccupancy / sum(virtualOccupancy); % normalize


    
    totalT = endT.(currEpoch) - startT.(currEpoch);
    nBins  = floor((totalT - binDur)/stepDur) + 1;

    binStarts               = startT.(currEpoch) + (0:nBins-1)*stepDur;
    binEnds                 = binStarts + binDur;
    binCenters.(currEpoch)  = binStarts + binDur/2;
    
    
    selectPBEs = cell(nBins, 1);
    nPBEs.(currEpoch) = zeros(nBins, 1);
    for ibin = 1:nBins
        selectPBEs{ibin} = find(peakT >= binStarts(ibin) & peakT < binEnds(ibin));
        nPBEs.(currEpoch)(ibin) = numel(selectPBEs{ibin});
    end
    

    % actual PBEs
    fprintf('\nActual PBEs')
    
    assemblyTuningCorrMat_time.(currEpoch)     = nan(nUnits, nUnits, nBins);
    assemblyTuningPFcorr_time.(currEpoch)      = nan(nUnits, nBins);
    assemblyTuningSpatialInfo_time.(currEpoch) = nan(nUnits, nBins);
    assemblyTuningPFKLdiv_time.(currEpoch)     = nan(nUnits, nBins);
    assemblyTuning_prefPos_time.(currEpoch)    = nan(nUnits, nBins);
    assemblyTuningPF_prefPosDist.(currEpoch)   = nan(nUnits, nBins);
    
    ifshuffle = 0;
    currAssemblyTunings = calculateAssemblyTuning_vs_time_v3(PBEInfo, spikes, selectPBEs, clusterQuality, ifshuffle);
    
    
    assemblyTunings_time.(currEpoch).data  = currAssemblyTunings;
    for ibin = 1:nBins
        
        assemblyTuningCorrMat_time.(currEpoch)(:, :, ibin) = corr(currAssemblyTunings(:, :, ibin)', spatialTunings_merge'); 
        assemblyTuningPFcorr_time.(currEpoch)(:, ibin)     = diag(assemblyTuningCorrMat_time.(currEpoch)(:, :, ibin)); % the correlation of PF and assembly tuning for a given unit

%         meanFR = sum(repmat(virtualOccupancy', [nUnits 1]).* currAssemblyTunings(:, :, ibin), 2) ./ sum(virtualOccupancy); % The calculation of mean firing rate could be a bit questionable
%         assemblyTuningSpatialInfo_time.(currEpoch)(:, ibin) = sum(repmat(virtualOccupancy', [nUnits 1]).* (currAssemblyTunings(:,:,ibin) ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings(:,:,ibin) ./ repmat(meanFR, [1 nPosBins])), 2);
    
%         KLdiv = calKLDivergence(currAssemblyTunings(:, :, ibin)', spatialTunings_merge');
%         assemblyTuningPFKLdiv_time.(currEpoch)(:, ibin) = diag(KLdiv);

%         [assemblyTuningPF_prefPosDist.(currEpoch)(:, ibin), assemblyTuning_prefPos_time.(currEpoch)(:, ibin)] = calDistPrefPosition(currAssemblyTunings(:, :, ibin), PF_prefPos);
    
    end
   
end


fileName = [sessionName '.assemblyTuning_vs_time_Jan2020.mat'];

save(fullfile(storagePath, fileName), 'nPBEs', 'binCenters', ...
        'assemblyTunings_time', 'assemblyTuningCorrMat_time', 'assemblyTuningPFcorr_time', ...
        'assemblyTuningSpatialInfo_time', ...
        'assemblyTuningPFKLdiv_time', ...
        'assemblyTuning_prefPos_time', 'assemblyTuningPF_prefPosDist', '-v7.3')



end

%% sub-functions


function KLdiv = calKLDivergence(learnedTunings, spatialTunings)


% the individual learned tunings and spatial tunings should be in columns
% this function works well when when the learnedTunings matrix contains multipels instance of learned tunings (for examples if they belong to different time windows or different instances of unit identity shuffle)


nPosBins = size(spatialTunings, 1);

spatialTunings = spatialTunings ./ repmat(sum(spatialTunings, 1), [nPosBins 1]);
spatialTunings = spatialTunings + eps;

learnedTunings = learnedTunings ./ repmat(sum(learnedTunings, 1), [nPosBins 1]);
learnedTunings = learnedTunings + eps;


nSTs = size(spatialTunings, 2);
nLTs = size(learnedTunings, 2); % number of learned tunings


KLdiv = nan(nLTs, nSTs);
for ii = 1:nSTs
    
    currSpatialTuning = spatialTunings(:, ii);
   
    spatialTuningTerm = repmat(currSpatialTuning, [1 nLTs]);
    KLdiv(:, ii) = sum(spatialTuningTerm .* (log(spatialTuningTerm) - log(learnedTunings)), 1);

end


end


function [prefPosDist, asTuning_prefPos] = calDistPrefPosition(assemblyTunings, PF_prefPos)  
    
nUnits     = numel(PF_prefPos);
nInstances = size(assemblyTunings, 3);
nPosBins   = size(assemblyTunings, 2);

prefPosDist      = nan(nUnits, nInstances);
asTuning_prefPos = nan(nUnits, nInstances);

for inst = 1:nInstances
   currAsTunings = assemblyTunings(:,:, inst);
   
   [~, asTuning_prefPos(:, inst)] = max(currAsTunings, [], 2);
   prefPosDist(:, inst)           = abs(asTuning_prefPos(:, inst) - PF_prefPos(:));
   prefPosDist(:, inst)           = prefPosDist(:, inst)/nPosBins;
   
end

end

