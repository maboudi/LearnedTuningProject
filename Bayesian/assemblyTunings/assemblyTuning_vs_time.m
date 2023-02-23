function assemblyTuning_vs_time(sessionNumber)


% we might need to consider only stable units when calculating the assembly
% tuning for a specific unit

sz = getenv('SLURM_CPUS_PER_TASK');
p  = parpool('local', str2double(sz));


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

    virtualOccupancy.(currPeriod) = mean(posteriorProbMatrix, 2);
    virtualOccupancy.(currPeriod) = virtualOccupancy.(currPeriod) / sum(virtualOccupancy.(currPeriod)); % normalize


    
    totalT = endT.(currPeriod) - startT.(currPeriod);
    nBins  = floor((totalT - binDur)/stepDur) + 1;

    binStarts   = startT.(currPeriod) + (0:nBins-1)*stepDur;
    binEnds     = binStarts + binDur;
    binCenters.(currPeriod) = binStarts + binDur/2;

    
    nPBEs.(currPeriod)                          = zeros(nBins, 1);
    
    assemblyTunings_time.(currPeriod)           = zeros(nUnits, nPosBins, nBins);        
    assemblyTuningCorrMat_time.(currPeriod)     = zeros(nUnits, nUnits, nBins);
    assemblyTuningPFcorr_time.(currPeriod)      = zeros(nUnits, nBins);
    assemblyTuningSpatialInfo_time.(currPeriod) = zeros(nUnits, nBins);   
    
    assemblyTunings_time_zscore.(currPeriod)           = zeros(nUnits, nPosBins, nBins); 
    assemblyTuningPFcorr_time_zscore.(currPeriod)      = zeros(nUnits, nBins);
    assemblyTuningSpatialInfo_time_zscore.(currPeriod) = zeros(nUnits, nBins);   

    
    for ibin = 1:nBins
       
        idx = peakT >= binStarts(ibin) & peakT < binEnds(ibin);
        
        
        if ~isempty(idx)
            currPBEs = PBEInfo(idx);
            nPBEs.(currPeriod)(ibin) = numel(currPBEs); 

            
            
            % actual PBEs
            fprintf('\nCalculating asTunings bin %d out of %d - %s..', ibin, nBins, currPeriod)

            fprintf('\nActual PBEs')
            selectPBEs = [];
            ifShuffle  = 0;

            currAssemblyTunings = calculateAssemblyTuningV6(currPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);


            assemblyTunings_time.(currPeriod)(:, :, ibin)       = currAssemblyTunings;

            assemblyTuningCorrMat_time.(currPeriod)(:, :, ibin) = corr(currAssemblyTunings', spatialTunings_merge'); 
            assemblyTuningPFcorr_time.(currPeriod)(:, ibin)     = diag(assemblyTuningCorrMat_time.(currPeriod)(:, :, ibin)); % the correlation of PF and assembly tuning for a given unit


            meanFR = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* currAssemblyTunings, 2) ./ sum(virtualOccupancy.(currPeriod)); % The calculation of mean firing rate could be a bit questionable
            assemblyTuningSpatialInfo_time.(currPeriod)(:, ibin) = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* (currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);



            % unit ID shuffle PBEs
            fprintf('\nUnit ID shuffle PBEs')

            ifShuffle = 1;
            currAssemblyTunings = zeros(nUnits, nPosBins, nUIshuffles); % we are not going to store this variable, only need it for z-scoring 
            currPFcorr          = zeros(nUnits, nUIshuffles);
            currSpatialInfo     = zeros(nUnits, nUIshuffles);

            currOccupancy       = virtualOccupancy.(currPeriod);


            for seg = 1:nSeg
                fprintf('\nSegment %d out of %d', seg, nseg)
                parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)

                    shuffleAssemblyTunings = calculateAssemblyTuningV6(currPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);
                    currAssemblyTunings(:, :, i_ui) = shuffleAssemblyTunings;

                    currCorrMat         = corr(shuffleAssemblyTunings', spatialTunings_merge'); 
                    currPFcorr(:, i_ui) = diag(currCorrMat);

                    meanFR                   = sum(repmat(currOccupancy', [nUnits 1]).* shuffleAssemblyTunings, 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
                    currSpatialInfo(:, i_ui) = sum(repmat(currOccupancy', [nUnits 1]).* (shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);

                end
            end

            % comparing actual assembly tunings with the assembly tunings based on
            % unit ID shuffles
            assemblyTunings_time_zscore.(currPeriod)(:, :, ibin)  = (assemblyTunings_time.(currPeriod)(:, :, ibin) - mean(currAssemblyTunings, 3)) ./ std(currAssemblyTunings, [], 3);

            assemblyTuningPFcorr_time_zscore.(currPeriod)(:, ibin)      = (assemblyTuningPFcorr_time.(currPeriod)(:, ibin) - mean(currPFcorr, 2)) ./ std(currPFcorr, [], 2);
            assemblyTuningSpatialInfo_time_zscore.(currPeriod)(:, ibin) = (assemblyTuningSpatialInfo_time.(currPeriod)(:, ibin) - mean(currSpatialInfo, 2)) ./ std(currSpatialInfo, [], 2);
        end
    end

end


fileName = [sessionName '.assemblyTuning_vs_time.mat'];
save(fullfile(storagePath, fileName), 'nPBEs',...
    'assemblyTunings_time', 'assemblyTuningCorrMat_time', 'assemblyTuningPFcorr_time', 'assemblyTuningSpatialInfo_time', ...
    'assemblyTunings_time_zscore', 'assemblyTuningPFcorr_time_zscore', 'assemblyTuningSpatialInfo_time_zscore', 'binCenters')



