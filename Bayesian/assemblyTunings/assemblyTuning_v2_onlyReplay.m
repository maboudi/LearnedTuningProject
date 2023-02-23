function assemblyTuning_v2_onlyReplay(sessionNumber)


% we might need to consider only stable units when calculating the assembly
% tuning for a specific unit

% sz = getenv('SLURM_CPUS_PER_TASK');
% 
% theCluster = parcluster('local');
% 
% JobFolder = sprintf('/home/kmaboudi/.matlab/trash/job_%s', sessionNumber);
% mkdir(JobFolder)
% theCluster.JobStorageLocation = JobFolder;
% 
% p = parpool(theCluster, str2double(sz));



% parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name


basePath    = fullfile(parentDir, sessionName);
storagePath = fullfile(basePath, 'assemblyTunings');
mkdir(storagePath)


load(fullfile(basePath, 'BayesianDecoding', [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')

clusterQuality = clusterQuality; 


% clusterQuality = [];
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


epochs = fileInfo.behavior.time;

% we are going to analyze assembly tuning separately for each epoch, so we
% need to reorganize the PBEs


% PBEInfo = PBEInfo(acceptedIdx); % only PBEs qulaified in terms of duration, number of participant units, the brain satte during which they occurred
PBEInfo = PBEInfo_replayScores;

epochNames = {'pre', 'run', 'post'};



for pbe = 1:numel(PBEInfo)
    PBEInfo(pbe).jd = PBEInfo(pbe).jumpDist.maxJump;
end
allJDs = [PBEInfo.jd];

JD_thresh = prctile(allJDs, 50);


for iEpoch = 1:numel(epochNames) 
    currEpoch = epochNames{iEpoch};
    
    PBEs.(currEpoch)  = PBEInfo(strcmp({PBEInfo.epoch}, currEpoch));
    
    peaks = [PBEs.(currEpoch).peakT];
    switch iEpoch
        case 1
            idx = peaks > (epochs(1,2)-4*60*60);
        case 2
            idx = true(numel(peaks), 1);
        case 3
            idx  = peaks < (epochs(3,1)+4*60*60);
    end
   
    PBEs.(currEpoch) = PBEs.(currEpoch)(idx);
    nPBEs.(currEpoch) = numel(PBEs.(currEpoch));
end



%% assembly Tunings considering all of the PBEs
% Spatial information and correlation with place fields of assesmbly tunings
% and corresponding z-scores by comparing the values againts distibutions calculated based on unit ID shuffle surrogates


% PF_peak_min = 1;
% minPBEsParticipation.pre  = 100;
% minPBEsParticipation.run  = 5;
% minPBEsParticipation.post = 100;
% 
% nUIshuffles = 100;
% nSeg        = 5;
% 
% nUIshuffles_perSeg = nUIshuffles/nSeg;

for iEpoch = 1:3
    
    
    currEpoch = epochNames{iEpoch};
    fprintf(['\nProcessing ' currEpoch ' ..'])
    
    
    currPBEs = PBEs.(currEpoch);
    
    
%     % calculating the participatino of the units during the PBEs
%     participation = nan(nUnits, nPBEs.(currEpoch));
%     for ipbe = 1:nPBEs.(currEpoch)
%         participation(:, ipbe) = sum(currPBEs(ipbe).fr_20msbin, 2) > 0;            
%     end
% 
%     nParticipatedPBEs.(currEpoch) = sum(participation, 2);
%     activeUnits.(currEpoch)       = find(max(spatialTunings_merge, [], 2) > PF_peak_min & nParticipatedPBEs.(currEpoch) > minPBEsParticipation.(currEpoch));

    
    
    
    posteriorProbMatrix = cell(nPBEs.(currEpoch), 1);
    for pbe = 1:nPBEs.(currEpoch) 
        posteriorProbMatrix{pbe} = currPBEs(pbe).posteriorProbMat;
    end
    posteriorProbMatrix = cell2mat(posteriorProbMatrix');
    
    virtualOccupancy.(currEpoch) = mean(posteriorProbMatrix, 2);
    virtualOccupancy.(currEpoch) = virtualOccupancy.(currEpoch) / sum(virtualOccupancy.(currEpoch)); % normalize
    


    
    % actual PBEs
%     fprintf('\nCalculating the assembly tunings based on actual PBEs ..')
%     
%     selectPBEs = [];
%     ifShuffle  = 0;
%     currAssemblyTunings = calculateAssemblyTuningV8(currPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);
%     
%     
%     assemblyTunings.(currEpoch).data = currAssemblyTunings;
%     
%     assemblyTuningCorrMat.(currEpoch).data = corr(currAssemblyTunings', spatialTunings_merge', 'type', 'pearson'); 
%     assemblyTuningPFcorr.(currEpoch).data  = diag(assemblyTuningCorrMat.(currEpoch).data); % the correlation of PF and assembly tuning for a given unit
%     
%     meanFR = sum(repmat(virtualOccupancy.(currEpoch)', [nUnits 1]).* currAssemblyTunings, 2) ./ sum(virtualOccupancy.(currEpoch)); % The calculation of mean firing rate could be a bit questionable
%     assemblyTuningSpatialInfo.(currEpoch).data = sum(repmat(virtualOccupancy.(currEpoch)', [nUnits 1]).* (currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);
% 
%     KLdiv = calKLDivergence(currAssemblyTunings', spatialTunings_merge');
%     assemblyTuningPFKLdiv.(currEpoch).data  = diag(KLdiv);
%     
%     
%     [assemblyTuningPF_prefPosDist.(currEpoch).data, assemblyTuning_prefPos.(currEpoch).data] = calDistPrefPosition(assemblyTunings.(currEpoch).data, PF_prefPos);
%     
%     
% 
%     % unit ID shuffle PBEs
%     fprintf('\nCalculating assembly tunings based on unit ID shuffle surrogate PBEs ..')
%     
%     currAssemblyTunings = nan(nUnits, nPosBins, nUIshuffles);
%     currCorrMat         = nan(nUnits, nUnits, nUIshuffles);
%     currPFcorr          = nan(nUnits, nUIshuffles);
%     currSpatialInfo     = nan(nUnits, nUIshuffles);
%     currKLdiv           = nan(nUnits, nUIshuffles);
%     currPrefPos         = nan(nUnits, nUIshuffles);
%     currPrefPosDist     = nan(nUnits, nUIshuffles);
%     
%     
%     currOccupancy       = virtualOccupancy.(currEpoch);
%     
%     selectPBEs          = []; % all PBEs are included
%     ifShuffle           = 1;
% 
%     
%     for seg = 1:nSeg
%         seg
%         parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)
% 
%             shuffleAssemblyTunings = calculateAssemblyTuningV8(currPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);
% 
%             currAssemblyTunings(:, :, i_ui) = shuffleAssemblyTunings;
% 
%             currCorrMat(:, :, i_ui)  = corr(shuffleAssemblyTunings', spatialTunings_merge', 'type', 'pearson'); 
%             currPFcorr(:, i_ui)      = diag(currCorrMat(:, :, i_ui));
% 
%             meanFR                   = sum(repmat(currOccupancy', [nUnits 1]).* shuffleAssemblyTunings, 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
%             currSpatialInfo(:, i_ui) = sum(repmat(currOccupancy', [nUnits 1]).* (shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);
%             
%             KLdiv                    = calKLDivergence(shuffleAssemblyTunings', spatialTunings_merge');
%             currKLdiv(:, i_ui)       = diag(KLdiv);
%             
%             [currPrefPosDist(:, i_ui), currPrefPos(:, i_ui)] = calDistPrefPosition(shuffleAssemblyTunings, PF_prefPos);
%             
%         end
%     end
%     
%     assemblyTunings.(currEpoch).ui               = currAssemblyTunings;
%     assemblyTuningCorrMat.(currEpoch).ui         = currCorrMat;
%     assemblyTuningPFcorr.(currEpoch).ui          = currPFcorr;
%     assemblyTuningSpatialInfo.(currEpoch).ui     = currSpatialInfo;
%     assemblyTuningPFKLdiv.(currEpoch).ui         = currKLdiv;
%     assemblyTuning_prefPos.(currEpoch).ui        = currPrefPos;
%     assemblyTuningPF_prefPosDist.(currEpoch).ui  = currPrefPosDist;
%     
%     
%     % comparing actual assembly tunings with the assembly tunings based on
%     % unit ID shuffles
%     
%     assemblyTunings_zscore.(currEpoch)              = calZscore(assemblyTunings.(currEpoch), 3);
%     assemblyTuningPFcorr_zscore.(currEpoch)         = calZscore(assemblyTuningPFcorr.(currEpoch), 2);
%     assemblyTuningSpatialInfo_zscore.(currEpoch)    = calZscore(assemblyTuningSpatialInfo.(currEpoch), 2);
%     assemblyTuningPFKLdiv_zscore.(currEpoch)        = calZscore(assemblyTuningPFKLdiv.(currEpoch), 2);
%     assemblyTuningPF_prefPosDist_zscore.(currEpoch) = calZscore(assemblyTuningPF_prefPosDist.(currEpoch), 2);
    
    
end


% fileName = [sessionName '.assemblyTunings_allPBEs_Dec2021.mat'];
% 
% save(fullfile(storagePath, fileName), ...
%     'nPBEs', ...
%     'nParticipatedPBEs', ...
%     'activeUnits', ...
%     'virtualOccupancy', ...
%     'assemblyTunings', ...
%     'assemblyTuningCorrMat', ...
%     'assemblyTuningPFcorr', ...
%     'assemblyTuningSpatialInfo', ...
%     'assemblyTuningPFKLdiv', ...
%     'assemblyTuning_prefPos', ...
%     'assemblyTuningPF_prefPosDist' , ...
%     'assemblyTunings_zscore', ...
%     'assemblyTuningPFcorr_zscore', ...
%     'assemblyTuningSpatialInfo_zscore', ...
%     'assemblyTuningPFKLdiv_zscore', ...
%     'assemblyTuningPF_prefPosDist_zscore')


% return



%% How assembly tunigns varies as a function of BD replay score


nPBEs_subset = struct('pre', [], 'run', [], 'post', []);

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
           

nUIshuffles = 100;
nSeg        = 5;

nUIshuffles_perSeg = nUIshuffles/nSeg;


for irsm = 5
    
    replayScoreMethod = replayScoreMethods{irsm};
    replayScoreMethod_name = replayScoreMethods_fullName{irsm};
    
    fprintf(['\nReplay score method: ' replayScoreMethod_name ' ..'])
    
    
    for iEpoch = 1:3 % behavioral epochs: pre, run, post
        
        currEpoch = epochNames{iEpoch};
        
        fprintf(['\nepoch ' currEpoch ' ..'])
        
        
        currBDscores = [PBEs.(currEpoch).(replayScoreMethod)];
        
        %%% added Aug 5th 2022
        currJumpDistances = [PBEs.(currEpoch).jd];




        

        nPBEs_subset.(currEpoch).(replayScoreMethod) = nan(4,1);


        assemblyTunings_sub.(currEpoch).(replayScoreMethod).data        = nan(nUnits, nPosBins, 4);        
        assemblyTunings_sub.(currEpoch).(replayScoreMethod).ui          = cell(4,1);

        asTuningSpatialInfo_sub.(currEpoch).(replayScoreMethod).data    = nan(nUnits, 4);
        asTuningPFcorr_sub.(currEpoch).(replayScoreMethod).data         = nan(nUnits, 4);
        asTuningCorrMat_sub.(currEpoch).(replayScoreMethod).data        = nan(nUnits, nUnits, 4); 
        asTuningPFKLdiv_sub.(currEpoch).(replayScoreMethod).data        = nan(nUnits, 4);
        asTuning_prefPos_sub.(currEpoch).(replayScoreMethod).data       = nan(nUnits, 4);
        asTuningPF_prefPosDist_sub.(currEpoch).(replayScoreMethod).data = nan(nUnits, 4);
        
        
        asTuningSpatialInfo_sub.(currEpoch).(replayScoreMethod).ui      = nan(nUnits, 4, nUIshuffles);
        asTuningPFcorr_sub.(currEpoch).(replayScoreMethod).ui           = nan(nUnits, 4, nUIshuffles);
        asTuningPFKLdiv_sub.(currEpoch).(replayScoreMethod).ui          = nan(nUnits, 4, nUIshuffles);
        asTuning_prefPos_sub.(currEpoch).(replayScoreMethod).ui         = nan(nUnits, 4, nUIshuffles);
        asTuningPF_prefPosDist_sub.(currEpoch).(replayScoreMethod).ui   = nan(nUnits, 4, nUIshuffles);
 


        for iPBEset = 1:4


            fprintf('\nCalculating assembly tunings based on PBEs in quartile %d of replay score ..', iPBEset)


%             selectPBEs = find(currBDscores' > (iPBEset-1)*25 & currBDscores' <= iPBEset*25); 

            selectPBEs = find(currBDscores' > (iPBEset-1)*25 & currBDscores' <= iPBEset*25 & currJumpDistances' < JD_thresh); 

            

            nPBEs_subset.(currEpoch).(replayScoreMethod)(iPBEset) = numel(selectPBEs);
            currbinnedPBEs = PBEs.(currEpoch);
            


            % calculate the assembly tunings based on the current PBEs
            
            ifShuffle = 0;
            currAssemblyTunings = calculateAssemblyTuningV7(currbinnedPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);

            assemblyTunings_sub.(currEpoch).(replayScoreMethod).data(:, :, iPBEset) = currAssemblyTunings;

%             meanFR  = sum(repmat(virtualOccupancy.(currEpoch)', [nUnits 1]).* currAssemblyTunings, 2) ./ sum(virtualOccupancy.(currEpoch));
%             asTuningSpatialInfo_sub.(currEpoch).(replayScoreMethod).data(:, iPBEset) = sum(repmat(virtualOccupancy.(currEpoch)', [nUnits 1]).* (currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);

            asTuningCorrMat_sub.(currEpoch).(replayScoreMethod).data(:, :, iPBEset) = corr(currAssemblyTunings', spatialTunings_merge', 'type', 'pearson');
            asTuningPFcorr_sub.(currEpoch).(replayScoreMethod).data(:, iPBEset)     = diag(asTuningCorrMat_sub.(currEpoch).(replayScoreMethod).data(:, :, iPBEset));
            
            
%             KLdiv = calKLDivergence(currAssemblyTunings', spatialTunings_merge');
%             asTuningPFKLdiv_sub.(currEpoch).(replayScoreMethod).data(:, iPBEset)    = diag(KLdiv);
%             
%             
%             [asTuningPF_prefPosDist_sub.(currEpoch).(replayScoreMethod).data(:, iPBEset), ...
%                 asTuning_prefPos_sub.(currEpoch).(replayScoreMethod).data(:, iPBEset)] = calDistPrefPosition(currAssemblyTunings, PF_prefPos);


% 
%             % unit ID shuffle PBEs
%             
%             ifShuffle = 1;
%             currAssemblyTunings = nan(nUnits, nPosBins, nUIshuffles); % we are not going to store this variable, only need it for z-scoring 
%             currPFcorr          = nan(nUnits, nUIshuffles);
%             currSpatialInfo     = nan(nUnits, nUIshuffles);
%             currKLdiv           = nan(nUnits, nUIshuffles);
%             currPrefPos         = nan(nUnits, nUIshuffles);
%             currPrefPosDist     = nan(nUnits, nUIshuffles);
%             
%             
%             currOccupancy       = virtualOccupancy.(currEpoch);
% 
% 
%             fprintf('\nCalculating assembly tunings based on unit ID shuffle surrogate PBEs ..')
% 
%             for seg = 1:nSeg
%                 seg
%                 parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)
% 
%                     shuffleAssemblyTunings          = calculateAssemblyTuningV8(currbinnedPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);
%                     currAssemblyTunings(:, :, i_ui) = shuffleAssemblyTunings;
% 
%                     currCorrMat         = corr(shuffleAssemblyTunings', spatialTunings_merge', 'type', 'pearson'); 
%                     currPFcorr(:, i_ui) = diag(currCorrMat);
% 
%                     meanFR                   = sum(repmat(currOccupancy', [nUnits 1]).* shuffleAssemblyTunings, 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
%                     currSpatialInfo(:, i_ui) = sum(repmat(currOccupancy', [nUnits 1]).* (shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);
%                     
%                     
%                     KLdiv                    = calKLDivergence(shuffleAssemblyTunings', spatialTunings_merge');
%                     currKLdiv(:, i_ui)       = diag(KLdiv);
%                     
%                     [currPrefPosDist(:, i_ui), currPrefPos(:, i_ui)] = calDistPrefPosition(shuffleAssemblyTunings, PF_prefPos);
% 
%                 end
%             end
% 
%             assemblyTunings_sub.(currEpoch).(replayScoreMethod).ui{iPBEset}              = currAssemblyTunings;
%             asTuningPFcorr_sub.(currEpoch).(replayScoreMethod).ui(:, iPBEset, :)         = permute(currPFcorr, [1 3 2]);
%             asTuningSpatialInfo_sub.(currEpoch).(replayScoreMethod).ui(:, iPBEset, :)    = permute(currSpatialInfo, [1 3 2]);
%             asTuningPFKLdiv_sub.(currEpoch).(replayScoreMethod).ui(:, iPBEset, :)        = permute(currKLdiv, [1 3 2]);
%             asTuning_prefPos_sub.(currEpoch).(replayScoreMethod).ui(:, iPBEset, :)       = permute(currPrefPos, [1 3 2]);
%             asTuningPF_prefPosDist_sub.(currEpoch).(replayScoreMethod).ui(:, iPBEset, :) = permute(currPrefPosDist, [1 3 2]);



            % comparing actual assembly tunings with the assembly tunings based on
            % unit ID shuffles
            
%             dim = 3;
%             asTuningPFcorr_sub_zscore.(currEpoch).(replayScoreMethod)         = calZscore(asTuningPFcorr_sub.(currEpoch).(replayScoreMethod), dim);
%             asTuningSpatialInfo_sub_zscore.(currEpoch).(replayScoreMethod)    = calZscore(asTuningSpatialInfo_sub.(currEpoch).(replayScoreMethod), dim);
%             asTuningPFKLdiv_sub_zscore.(currEpoch).(replayScoreMethod)        = calZscore(asTuningPFKLdiv_sub.(currEpoch).(replayScoreMethod), dim);
%             asTuningPF_prefPosDist_sub_zscore.(currEpoch).(replayScoreMethod) = calZscore(asTuningPF_prefPosDist_sub.(currEpoch).(replayScoreMethod), dim);
%                  

        end

    end
end


fileName = [sessionName '.assemblyTuning_vs_replayScores_Lthresh1e_3_wcts_lowMaxJump.mat'];

save(fullfile(storagePath, fileName), ...
    'nPBEs_subset', ...
    'assemblyTunings_sub', ...
    'asTuningCorrMat_sub', ...
    'asTuningPFcorr_sub', ...
    'asTuningSpatialInfo_sub', ...
    'asTuningPFKLdiv_sub', ...
    'asTuning_prefPos_sub', ...
    'asTuningPF_prefPosDist_sub')%, ...
%     'asTuningPFcorr_sub_zscore', ...
%     'asTuningSpatialInfo_sub_zscore', ...
%     'asTuningPFKLdiv_sub_zscore', ...
%     'asTuningPF_prefPosDist_sub_zscore')


end

%% functions


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


function zscored = calZscore(variable, dim)

zscored  = (variable.data - nanmean(variable.ui, dim)) ./ nanstd(variable.ui, [], dim);

end
