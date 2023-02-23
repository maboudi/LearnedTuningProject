function assemblyTuning(sessionNumber)


sz = getenv('SLURM_CPUS_PER_TASK');
p  = parpool('local', str2num(sz));



addpath(genpath('/nfs/turbo/umms-kdiba/Kourosh/ReplayPreplayAnalyses'))

% cd('/home/kouroshmaboudi/Documents/HMM_project/assemblyTunings_greatLakes_reclusteredGrosmark_Hiro')
% currDir = '/home/kouroshmaboudi/Documents/HMM_project/assemblyTunings_greatLakes_reclusteredGrosmark_Hiro/';

currDir = '/nfs/turbo/umms-kdiba/Kourosh/assemblyTunings_greatLakes_reclusteredGrosmark_Hiro';
cd(currDir)


rr = dir;


% for sessionNumber = 1:8

% sessionNumber = 1;
    
sessionName = rr(sessionNumber+2).name;

% currSessionFolder = fullfile(currDir, sessionName, 'BayesianDecoding');
currSessionFolder = fullfile(currDir, sessionName);


cd(currSessionFolder)

% ll = dir;
% try
%     load(ll(4).name)
% catch
%     load(ll(3).name)
% end
load('allVariables.mat')
load([sessionName '.clusterQuality.mat']);

cq = clusterQuality;

unitIDs = [shanks' clus'];


if exist('spatialTunings', 'var')
    spatialTunings_merge = spatialTunings;
    spatialTunings_A     = spatialTunings;
    spatialTunings_B     = [];
elseif exist('spatialTunings_RL', 'var')
    
%     spatialTunings_merge = (spatialTunings_RL + spatialTunings_LR)/2;
    spatialTunings_merge = spatialTunings_biDir;
    spatialTunings_A     = spatialTunings_RL;
    spatialTunings_B     = spatialTunings_LR;
end


[nUnits, nPosBins] = size(spatialTunings_merge);


activeUnits  = find(sum(spatialTunings_merge, 2));
nActiveUnits = numel(activeUnits); 


[~, peakPos] = max(spatialTunings_merge(activeUnits, :), [], 2);
[~, sortInd] = sort(peakPos, 'descend');

sortedActiveUnitIndices = activeUnits(sortInd);




%% selected PBEs

binnedPBEs.PRE  = PREbinnedPBEs;
binnedPBEs.RUN  = RUNbinnedPBEs;
binnedPBEs.POST = POSTbinnedPBEs;


periodNames = {'PRE', 'RUN', 'POST'};

PBELen = struct('PRE', [], 'RUN', [], 'POST', []);

for iperiod = 1:3

    currentPBEs = binnedPBEs.(periodNames{iperiod});
    nPBEs       = size(currentPBEs, 1);

    PBELen.(periodNames{iperiod})       = zeros(nPBEs, 1);
    nFiringUnits.(periodNames{iperiod}) = zeros(nPBEs, 1);

    for pbe = 1:nPBEs

        PBELen.(periodNames{iperiod})(pbe)       = size(currentPBEs{pbe, 2}, 2);
        nFiringUnits.(periodNames{iperiod})(pbe) = numel(find(sum(currentPBEs{pbe, 2}, 2)));

    end
end

clear nPBEs



if isfield(behavior, 'time') % one of Hiro's sessions
    
   secondaryPBEs2.PRE  = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(1, 1) & secondaryPBEs(:, 2) < behavior.time(1, 2), :);
   secondaryPBEs2.RUN  = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(2, 1) & secondaryPBEs(:, 2) < behavior.time(2, 2), :);
   secondaryPBEs2.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(3, 1) & secondaryPBEs(:, 2) < behavior.time(3, 2), :); 
    
else % one of Grosmark's sessions
    
   secondaryPBEs2.PRE  = secondaryPBEs(secondaryPBEs(:, 1) > behavior.PREEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.PREEpoch(1, 2), :);
   secondaryPBEs2.RUN  = secondaryPBEs(secondaryPBEs(:, 1) > behavior.MazeEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.MazeEpoch(1, 2), :);
   secondaryPBEs2.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.POSTEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.POSTEpoch(1, 2), :);
   
end


acceptedEvts.PRE  = find(secondaryPBEs2.PRE(:, 6) == 1 & PBELen.PRE >= 4 & nFiringUnits.PRE >= 5); % 5:RippleOverlap, 6:NREM, 7:REM, 8:QW, 9:wake
acceptedEvts.RUN  = find(secondaryPBEs2.RUN(:, 5) == 1 & PBELen.RUN >= 4 & nFiringUnits.RUN >= 5);
acceptedEvts.POST = find(secondaryPBEs2.POST(:, 6) == 1 & PBELen.POST >= 4 & nFiringUnits.POST >= 5);


nPBEs.PRE  = numel(acceptedEvts.PRE);
nPBEs.RUN  = numel(acceptedEvts.RUN);
nPBEs.POST = numel(acceptedEvts.POST);



%% assembly Tunings considering all of the accepted PBEs
%% spatial information and correlation with place fields of assesmbly tunings
%% and their z-scores againts distibutions calculated based on unit ID shuffle surrogates


[assemblyTunings, assemblyTuningCorrMat, assemblyTuningPFcorr, assemblyTuningSpatialInfo, ...
    assemblyTunings_zscore, assemblyTuningPFcorr_zscore, assemblyTuningSpatialInfo_zscore, virtualOccupancy] = deal(struct('PRE', [], 'RUN', [], 'POST', []));

nUnitIDs        = 100;
nSeg            = 4;

nUnitIDs_perSeg = nUnitIDs/nSeg;

for iperiod = 1:3
    
    iperiod
    currPeriod = periodNames{iperiod};
   
    currPBEs   = binnedPBEs.(currPeriod);
   
    
    virtualOccupancy.(currPeriod) = mean(cell2mat(posteriorProbMatrix.(currPeriod).data'), 2);
    virtualOccupancy.(currPeriod) = virtualOccupancy.(currPeriod) / sum(virtualOccupancy.(currPeriod));

    
    fprintf('\nCalculating the assembly tunings based on actual PBEs ..')
    
    
    % actual PBEs
    currAssemblyTunings = calculateAssemblyTuningV5(currPBEs, spatialTunings_A, spatialTunings_B, acceptedEvts.(currPeriod), unitIDs, cq, 0);
%         currAssemblyTunings = calculateAssemblyTuningV3(currPBEs, spatialTunings_A, spatialTunings_B, acceptedEvts.(currPeriod), 0);
%                 currAssemblyTunings = calculateAssemblyTuningV4(currPBEs, spatialTunings_A, spatialTunings_B, acceptedEvts.(currPeriod), shanks',  0);


    assemblyTunings.(currPeriod).data = currAssemblyTunings(sortedActiveUnitIndices, :);
    
    assemblyTuningCorrMat.(currPeriod).data = corr(currAssemblyTunings(sortedActiveUnitIndices, :)', spatialTunings_merge(sortedActiveUnitIndices, :)'); 
    assemblyTuningPFcorr.(currPeriod).data  = diag(assemblyTuningCorrMat.(currPeriod).data);
    
    meanFR = sum(repmat(virtualOccupancy.(currPeriod)', [nActiveUnits 1]).* currAssemblyTunings(sortedActiveUnitIndices, :), 2) ./ sum(virtualOccupancy.(currPeriod)); % The calculation of mean firing rate could be a bit questionable
    assemblyTuningSpatialInfo.(currPeriod).data = sum(repmat(virtualOccupancy.(currPeriod)', [nActiveUnits 1]).* (currAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])), 2);

    
    
    % unit ID shuffle PBEs
    currAssemblyTunings = zeros(nActiveUnits, nPosBins, nUnitIDs);
    currCorrMat         = zeros(nActiveUnits, nActiveUnits, nUnitIDs);
    currPFcorr          = zeros(nActiveUnits, nUnitIDs);
    currSpatialInfo     = zeros(nActiveUnits, nUnitIDs);
    
    currOccupancy       = virtualOccupancy.(currPeriod);
    
    selectPBEs          = acceptedEvts.(currPeriod);
    sortedPFs           = spatialTunings_merge(sortedActiveUnitIndices, :);
    
    
    
    fprintf('\nCalculating the assembly tunings based on unit ID shuffle PBEs ..')
    
    for seg = 1:nSeg
        seg
        parfor i_ui = ((seg-1)*nUnitIDs_perSeg +1):(seg*nUnitIDs_perSeg)

            shuffleAssemblyTunings = calculateAssemblyTuningV5(currPBEs, spatialTunings_A, spatialTunings_B, selectPBEs, unitIDs, cq, 1);
            currAssemblyTunings(:, :, i_ui) = shuffleAssemblyTunings(sortedActiveUnitIndices, :);

            currCorrMat(:, :, i_ui)  = corr(shuffleAssemblyTunings(sortedActiveUnitIndices, :)', sortedPFs'); 
            currPFcorr(:, i_ui)      = diag(currCorrMat(:, :, i_ui));

            meanFR                   = sum(repmat(currOccupancy', [nActiveUnits 1]).* shuffleAssemblyTunings(sortedActiveUnitIndices, :), 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
            currSpatialInfo(:, i_ui) = sum(repmat(currOccupancy', [nActiveUnits 1]).* (shuffleAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])), 2);

        end
    end
    
    assemblyTunings.(currPeriod).ui           = currAssemblyTunings;
    assemblyTuningCorrMat.(currPeriod).ui     = currCorrMat;
    assemblyTuningPFcorr.(currPeriod).ui      = currPFcorr;
    assemblyTuningSpatialInfo.(currPeriod).ui = currSpatialInfo;
    
    
    
    % comparing actual assembly tunings with the assembly tunings based on
    % unit ID shuffles
    
    assemblyTunings_zscore.(currPeriod)           = (assemblyTunings.(currPeriod).data - mean(assemblyTunings.(currPeriod).ui, 3)) ./ std(assemblyTunings.(currPeriod).ui, [], 3);
    
    assemblyTuningPFcorr_zscore.(currPeriod)      = (assemblyTuningPFcorr.(currPeriod).data - mean(assemblyTuningPFcorr.(currPeriod).ui, 2)) ./ std(assemblyTuningPFcorr.(currPeriod).ui, [], 2);
    assemblyTuningSpatialInfo_zscore.(currPeriod) = (assemblyTuningSpatialInfo.(currPeriod).data - mean(assemblyTuningPFcorr.(currPeriod).ui, 2)) ./ std(assemblyTuningSpatialInfo.(currPeriod).ui, [], 2);
    
    
end

fileName = 'assemblyTunings_allPBEs.mat';
save(fullfile(currSessionFolder, fileName), 'nPBEs', 'acceptedEvts', 'virtualOccupancy', 'sortedActiveUnitIndices', ...
    'assemblyTunings', 'assemblyTuningCorrMat', 'assemblyTuningPFcorr', 'assemblyTuningSpatialInfo', ...
    'assemblyTunings_zscore', 'assemblyTuningPFcorr_zscore', 'assemblyTuningSpatialInfo_zscore')




%% How assembly tunigns varies as a function of BD replay score

[assemblyTunings_sub, asTuningCorrMat_sub, asTuningPFcorr_sub, asTuningSpatialInfo_sub, ...
    assemblyTunings_sub_zscore, asTuningPFcorr_sub_zscore, asTuningSpatialInfo_sub_zscore] = deal(struct('PRE', [], 'RUN', [], 'POST', [])); % assembly tunings calculated based on PBEs with replay scores in different quartiles of the replay score distributions

nPBEs_subset = struct('PRE', [], 'RUN', [], 'POST', []);


shuffleMethods       = {'wPBEtimeswap'; 'unitIDshuffle'};
replayScoringMethods = {'weightedCorr'; 'replayScore'};

nUnitIDs        = 100;
nSeg            = 4;

nUnitIDs_perSeg = nUnitIDs/nSeg;

for iseqmetric = 1:2 % radon transform or weighted correlation
    
    currReplayMethod  = replayScoringMethods{iseqmetric};
    
    
    for ishuffle = 1:2 % time swap or unit ID shuffle

        currShuffleMethod = shuffleMethods{ishuffle};
        
        
        for iperiod = 1:3 % behavioral epochs: PRE, RUN, POST

            
            currPeriod = periodNames{iperiod};
            
            currBDscores = BDseqscore.(currPeriod).data.(currShuffleMethod).(currReplayMethod).prctilescore;
            
            
            nPBEs_subset.(currPeriod).(currShuffleMethod).(currReplayMethod) = zeros(4,1);
            
            
            assemblyTunings_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)      = zeros(nActiveUnits, nPosBins, 4);        
            assemblyTunings_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)        = cell(4,1);
            
            asTuningSpatialInfo_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)  = zeros(nActiveUnits, 4);
            asTuningPFcorr_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)       = zeros(nActiveUnits, 4);
            asTuningCorrMat_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)      = zeros(nActiveUnits, nActiveUnits, 4);     
            
            
            asTuningSpatialInfo_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)    = zeros(nActiveUnits, 4, nUnitIDs);
            asTuningPFcorr_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)         = zeros(nActiveUnits, 4, nUnitIDs);

            
            virtualOccupancy.(currPeriod) = mean(cell2mat(posteriorProbMatrix.(currPeriod).data'), 2);
            virtualOccupancy.(currPeriod) = virtualOccupancy.(currPeriod) / sum(virtualOccupancy.(currPeriod));
          
            
            
            for iPBEset = 1:4

                
                fprintf('\nCalculating the assembly tunings based on PBEs in quartile %d of replay score ..', iPBEset)
                
                
                currPBEs = find(currBDscores' > (iPBEset-1)*25 & currBDscores' <= iPBEset*25 & ismember(1:numel(currBDscores), acceptedEvts.(currPeriod))); 
                
                nPBEs_subset.(currPeriod).(currShuffleMethod).(currReplayMethod)(iPBEset) = numel(currPBEs);
                currbinnedPBEs  = binnedPBEs.(currPeriod);

                

                
                % calculate the assembly tunings based on the current PBEs
                
                currAssemblyTunings = calculateAssemblyTuningV5(currbinnedPBEs, spatialTunings_A, spatialTunings_B, currPBEs, unitIDs, cq, 0);
                assemblyTunings_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, :, iPBEset) = currAssemblyTunings(sortedActiveUnitIndices, :);
                
                
                
                meanFR  = sum(repmat(virtualOccupancy.(currPeriod)', [nActiveUnits 1]).* currAssemblyTunings(sortedActiveUnitIndices, :), 2) ./ sum(virtualOccupancy.(currPeriod));
                asTuningSpatialInfo_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, iPBEset) = sum(repmat(virtualOccupancy.(currPeriod)', [nActiveUnits 1]).* (currAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])), 2);

                asTuningCorrMat_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, :, iPBEset) = corr(permute(currAssemblyTunings(sortedActiveUnitIndices, :), [2 1]), spatialTunings_merge(sortedActiveUnitIndices, :)');
                asTuningPFcorr_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, iPBEset)     = diag(asTuningCorrMat_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, :, iPBEset));
                
                
            

                
                % unit ID shuffle PBEs
                
                currAssemblyTunings = zeros(nActiveUnits, nPosBins, nUnitIDs); % we are not going to store this variable, only need it for z-scoring 
%                 currCorrMat         = zeros(nActiveUnits, nActiveUnits, nUnitIDs);
                currPFcorr          = zeros(nActiveUnits, nUnitIDs);
                currSpatialInfo     = zeros(nActiveUnits, nUnitIDs);
                
                currOccupancy       = virtualOccupancy.(currPeriod);

                sortedPFs           = spatialTunings_merge(sortedActiveUnitIndices, :);

                fprintf('\nCalculating the assembly tunings based on unit ID shuffle PBEs ..')
                
                
                for seg = 1:nSeg
                    seg
                    parfor i_ui = ((seg-1)*nUnitIDs_perSeg +1):(seg*nUnitIDs_perSeg)


                        shuffleAssemblyTunings = calculateAssemblyTuningV5(currbinnedPBEs, spatialTunings_A, spatialTunings_B, currPBEs, unitIDs, cq, 1);
                        currAssemblyTunings(:, :, i_ui) = shuffleAssemblyTunings(sortedActiveUnitIndices, :);

                        currCorrMat         = corr(shuffleAssemblyTunings(sortedActiveUnitIndices, :)', sortedPFs'); 
                        currPFcorr(:, i_ui) = diag(currCorrMat);

                        meanFR                   = sum(repmat(currOccupancy', [nActiveUnits 1]).* shuffleAssemblyTunings(sortedActiveUnitIndices, :), 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
                        currSpatialInfo(:, i_ui) = sum(repmat(currOccupancy', [nActiveUnits 1]).* (shuffleAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings(sortedActiveUnitIndices, :) ./ repmat(meanFR, [1 nPosBins])), 2);

                    end
                end
                
                assemblyTunings_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod){iPBEset}         = currAssemblyTunings;
                asTuningPFcorr_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)(:, iPBEset, :)      = permute(currPFcorr, [1 3 2]);
                asTuningSpatialInfo_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)(:, iPBEset, :) = permute(currSpatialInfo, [1 3 2]);
                
                
                
                % comparing actual assembly tunings with the assembly tunings based on
                % unit ID shuffles
                assemblyTunings_sub_zscore.(currPeriod).(currShuffleMethod).(currReplayMethod)(:, :, iPBEset) = (assemblyTunings_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, :, iPBEset) - mean(currAssemblyTunings, 3)) ./ std(currAssemblyTunings, [], 3);
                
                asTuningPFcorr_sub_zscore.(currPeriod).(currShuffleMethod).(currReplayMethod)(:, iPBEset)      = (asTuningPFcorr_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, iPBEset) - mean(asTuningPFcorr_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)(:, iPBEset, :), 3)) ./ std(asTuningPFcorr_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)(:, iPBEset, :), [], 3);
                asTuningSpatialInfo_sub_zscore.(currPeriod).(currShuffleMethod).(currReplayMethod)(:, iPBEset) = (asTuningSpatialInfo_sub.(currPeriod).data.(currShuffleMethod).(currReplayMethod)(:, iPBEset) - mean(asTuningSpatialInfo_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)(:, iPBEset, :), 3)) ./ std(asTuningSpatialInfo_sub.(currPeriod).ui.(currShuffleMethod).(currReplayMethod)(:, iPBEset, :), [], 3);


            end

        end
    end
end

fileName = 'assemblyTuning_vs_replayScores.mat';
save(fullfile(currSessionFolder, fileName), 'nPBEs_subset', ...
    'assemblyTunings_sub', 'asTuningCorrMat_sub', 'asTuningPFcorr_sub', 'asTuningSpatialInfo_sub', ...
    'assemblyTunings_sub_zscore', 'asTuningPFcorr_sub_zscore', 'asTuningSpatialInfo_sub_zscore')


end
