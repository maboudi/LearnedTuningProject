function assemblyTuning_smoothedSpikes_sweep(sessionNumber)


% we might need to consider only stable units when calculating the assembly
% tuning for a specific unit

% sz = getenv('SLURM_CPUS_PER_TASK');
% p  = parpool('local', str2double(sz));


% parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/concat_GreatLakes_datasets_temp';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name;


basePath    = fullfile(parentDir, sessionName);
storagePath = fullfile(basePath, 'assemblyTunings');
mkdir(storagePath)


load(fullfile(basePath, 'BayesianDecoding', [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')


spikes   = spikes_pyr;
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


epochs = fileInfo.behavior.time;


% we are going to analyze assembly tuning separately for each epoch, so we
% need to reorganize the PBEs


% PBEInfo = PBEInfo(acceptedIdx); % only PBEs qulaified in terms of duration, number of participant units, the brain satte during which they occurred

PBEInfo = PBEInfo_replayScores;

epochNames = {'pre', 'run', 'post'};


Qs  = 0:5;
nQs = numel(Qs); 


for iepoch = 1:3
    
    currEpoch = epochNames{iepoch};
    
    assemblyTunings.(currEpoch)       = zeros(nUnits, nPosBins, nQs);
    assemblyTuningCorrMat.(currEpoch) = zeros(nUnits, nUnits, nQs);
    assemblyTuningPFcorr.(currEpoch)  = zeros(nUnits, nQs);
    
    asTuning_PF_prefPosDist.(currEpoch)       = zeros(nUnits, nQs);
    asTuning_prefPos.(currEpoch)              = zeros(nUnits, nQs);
    asTuning_PF_prefPos_spearCorr.(currEpoch) = zeros(nQs, 1);
    asTuning_PF_prefPos_pval.(currEpoch)      = zeros(nQs, 1);
    
end


PF_peak_min = 1;
minPBEsParticipation.pre  = 100;
minPBEsParticipation.run  = 5;
minPBEsParticipation.post = 100;

for ismooth = 1:nQs

    
    if ismooth == 1 % no smoothing
        
        fprintf('\nCalculating the assembly tunings without smoothing the decoding spikes')
        for pbe = 1:numel(PBEInfo)
            PBEInfo(pbe).fr_20msbin_smoothed = PBEInfo(pbe).fr_20msbin;
        end

    else
        
        curr_Qs = Qs(ismooth);
        fprintf('\nCalculating the assembly tunings for decoding spikes smoothed by Qs = %d ..', curr_Qs)
        
        halfWidth = curr_Qs*4;
        gw = gausswindow(curr_Qs, halfWidth);


        for pbe = 1:numel(PBEInfo)

            s = PBEInfo(pbe).fr_20msbin;
            for iUnit = 1:nUnits
                s(iUnit, :) = conv(s(iUnit, :), gw, 'same');
            end

            PBEInfo(pbe).fr_20msbin_smoothed = s;
        end
        
        
    end


    for iepoch = 1:numel(epochNames) 
        currEpoch = epochNames{iepoch};

        PBEs.(currEpoch)  = PBEInfo(strcmp({PBEInfo.epoch}, currEpoch));
        if iepoch == 3
           peaks = [PBEs.(currEpoch).peakT];
            idx  = peaks < (epochs(3,1)+4*60*60);
            PBEs.(currEpoch) = PBEs.(currEpoch)(idx);
        end
        nPBEs.(currEpoch) = numel(PBEs.(currEpoch));
    end


    for iepoch = 1:3

        currEpoch = epochNames{iepoch};

        fprintf(['\n\nProcessing ' currEpoch ' ..'])

        currPBEs = PBEs.(currEpoch);


        selectPBEs = [];
        ifShuffle  = 0;
        currAssemblyTunings = calculateAssemblyTuningV6(currPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);



        assemblyTunings.(currEpoch)(:, :, ismooth) = currAssemblyTunings;

        assemblyTuningCorrMat.(currEpoch)(:, :, ismooth) = corr(currAssemblyTunings', spatialTunings_merge'); 
        assemblyTuningPFcorr.(currEpoch)(:, ismooth)  = diag(assemblyTuningCorrMat.(currEpoch)(:, :, ismooth)); % the correlation of PF and assembly tuning for a given unit\



        participation = zeros(nUnits, nPBEs.(currEpoch));
        for ipbe = 1:nPBEs.(currEpoch)
            participation(:, ipbe) = sum(currPBEs(ipbe).fr_20msbin, 2) > 0;            
        end

        nParticipatedPBEs.(currEpoch) = sum(participation, 2);
        activeUnits.(currEpoch) = find(max(spatialTunings_merge, [], 2) > PF_peak_min & nParticipatedPBEs.(currEpoch) > minPBEsParticipation.(currEpoch));

        [asTuning_PF_prefPosDist.(currEpoch)(:, ismooth), asTuning_prefPos.(currEpoch)(:, ismooth)] = calDistPrefPosition(assemblyTunings.(currEpoch)(:, :, ismooth), PF_prefPos, activeUnits.(currEpoch));
        [asTuning_PF_prefPos_spearCorr.(currEpoch)(ismooth), asTuning_PF_prefPos_pval.(currEpoch)(ismooth)] = calcorr(asTuning_prefPos.(currEpoch)(:, ismooth), PF_prefPos, activeUnits.(currEpoch));


    end

end


fileName = [sessionName '.assemblyTunings_allPBEs_smoothedPBEs.mat'];
save(fullfile(storagePath, fileName), 'nPBEs', 'activeUnits', 'assemblyTunings', 'assemblyTuningPFcorr', ...
    'asTuning_prefPos', 'asTuning_PF_prefPosDist', 'asTuning_PF_prefPos_spearCorr', 'asTuning_PF_prefPos_pval')

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

function [corrCoeff, pval] = calcorr(asTuning_prefPos, PF_prefPos, activeUnits)

nSets = size(asTuning_prefPos, 2);

corrCoeff = nan(nSets, 1);
pval      = nan(nSets, 1);

for n = 1:nSets   
    [corrCoeff(n), pval(n)] = corr(asTuning_prefPos(activeUnits, n), PF_prefPos(activeUnits), 'type', 'spearman');
end

end
