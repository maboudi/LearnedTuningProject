function assemblyTuning_v3(sessionNumber)


% for 10ms time windows instead of the usual 20ms


parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
% parentDir = '/data/concat_GreatLakes_datasets_temp';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name


sz = getenv('SLURM_CPUS_PER_TASK');

theCluster = parcluster('local');

JobFolder = sprintf('/home/kmaboudi/.matlab/trash2/job_%s', sessionNumber);
mkdir(JobFolder)
theCluster.JobStorageLocation = JobFolder;


p = parpool(theCluster, str2double(sz));



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
spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
end

epochs = fileInfo.behavior.time;

% we are going to analyze assembly tuning separately for each epoch, so we
% need to reorganize the PBEs


% PBEInfo = PBEInfo(acceptedIdx); % only PBEs qulaified in terms of duration, number of participant units, the brain satte during which they occurred
PBEInfo = PBEInfo_replayScores;


%% add the firing  rates in 50 ms time bins to the PBEInfo structure

for ipbe = 1:numel(PBEInfo)
    
    fr_1msbin = PBEInfo(ipbe).fr_1msbin;
    
    firstBin = find(sum(fr_1msbin, 1) > 0 , 1, 'first');
    lastBin  = find(sum(fr_1msbin, 1) > 0 , 1, 'last');
    
    fr_1msbin = fr_1msbin(:, firstBin:lastBin);
    
    nTimeBin = round(size(fr_1msbin, 2)/10);
    
    fr_10msbin = zeros(nUnits, nTimeBin);
    for ibin = 1:nTimeBin
        fr_10msbin(:, ibin) = sum(fr_1msbin(:, (ibin-1)*50+1 : min(ibin*50, size(fr_1msbin, 2))), 2);
    end
    
    PBEInfo(ipbe).fr_10msbin = fr_10msbin;
    
end



periodNames = {'pre', 'run', 'post'};

for iperiod = 1:numel(periodNames) 
    currPeriod = periodNames{iperiod};
    
    PBEs.(currPeriod)  = PBEInfo(strcmp({PBEInfo.epoch}, currPeriod));
    if iperiod == 3
       peaks = [PBEs.(currPeriod).peakT];
        idx  = peaks < (epochs(3,1)+4*60*60);
        PBEs.(currPeriod) = PBEs.(currPeriod)(idx);
    end
    nPBEs.(currPeriod) = numel(PBEs.(currPeriod));
end



%% assembly Tunings considering all of the PBEs
% Spatial information and correlation with place fields of assesmbly tunings
% and corresponding z-scores by comparing the values againts distibutions calculated based on unit ID shuffle surrogates


nUIshuffles = 100;
nSeg        = 5;

nUIshuffles_perSeg = nUIshuffles/nSeg;

for iperiod = 1:3
    
    
    currPeriod = periodNames{iperiod};
    
    fprintf(['/nProcessing ' currPeriod ' ..'])
    
    
    currPBEs = PBEs.(currPeriod);
    
    
    posteriorProbMatrix = cell(nPBEs.(currPeriod), 1);
    for pbe = 1:nPBEs.(currPeriod) 
        posteriorProbMatrix{pbe} = currPBEs(pbe).posteriorProbMat;
    end
    posteriorProbMatrix = cell2mat(posteriorProbMatrix');
    
    
    virtualOccupancy.(currPeriod) = mean(posteriorProbMatrix, 2);
    virtualOccupancy.(currPeriod) = virtualOccupancy.(currPeriod) / sum(virtualOccupancy.(currPeriod)); % normalize
    
    
%     
%     
%     % actual PBEs
%     
%     fprintf('\nCalculating the assembly tunings based on actual PBEs ..')
%     
%     selectPBEs = [];
%     ifShuffle  = 0;
%     currAssemblyTunings = calculateAssemblyTuningV7(currPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);
%     
%     
%     assemblyTunings.(currPeriod).data = currAssemblyTunings;
%     
%     assemblyTuningCorrMat.(currPeriod).data = corr(currAssemblyTunings', spatialTunings_merge'); 
%     assemblyTuningPFcorr.(currPeriod).data  = diag(assemblyTuningCorrMat.(currPeriod).data); % the correlation of PF and assembly tuning for a given unit
%     
%     
%     meanFR = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* currAssemblyTunings, 2) ./ sum(virtualOccupancy.(currPeriod)); % The calculation of mean firing rate could be a bit questionable
%     assemblyTuningSpatialInfo.(currPeriod).data = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* (currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);

    
    
    
    % unit ID shuffle PBEs
%     
%     fprintf('\nCalculating assembly tunings based on unit ID shuffle surrogate PBEs ..')
%     
%     currAssemblyTunings = zeros(nUnits, nPosBins, nUIshuffles);
%     currCorrMat         = zeros(nUnits, nUnits, nUIshuffles);
%     currPFcorr          = zeros(nUnits, nUIshuffles);
%     currSpatialInfo     = zeros(nUnits, nUIshuffles);
%     
%     currOccupancy       = virtualOccupancy.(currPeriod);
%     
%     selectPBEs          = []; % all PBEs are included
%     ifShuffle           = 1;
% 
%     for seg = 1:nSeg
%         seg
%         parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)
%             i_ui
% 
%             shuffleAssemblyTunings = calculateAssemblyTuningV7(currPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);
% 
%             currAssemblyTunings(:, :, i_ui) = shuffleAssemblyTunings;
% 
%             currCorrMat(:, :, i_ui)  = corr(shuffleAssemblyTunings', spatialTunings_merge'); 
%             currPFcorr(:, i_ui)      = diag(currCorrMat(:, :, i_ui));
% 
%             meanFR                   = sum(repmat(currOccupancy', [nUnits 1]).* shuffleAssemblyTunings, 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
%             currSpatialInfo(:, i_ui) = sum(repmat(currOccupancy', [nUnits 1]).* (shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);
% 
%         end
%     end
%     
%     assemblyTunings.(currPeriod).ui           = currAssemblyTunings;
%     assemblyTuningCorrMat.(currPeriod).ui     = currCorrMat;
%     assemblyTuningPFcorr.(currPeriod).ui      = currPFcorr;
%     assemblyTuningSpatialInfo.(currPeriod).ui = currSpatialInfo;
%     
%     
%     
%     % comparing actual assembly tunings with the assembly tunings based on
%     % unit ID shuffles
%     
%     assemblyTunings_zscore.(currPeriod)           = (assemblyTunings.(currPeriod).data - mean(assemblyTunings.(currPeriod).ui, 3)) ./ std(assemblyTunings.(currPeriod).ui, [], 3);
%     assemblyTuningPFcorr_zscore.(currPeriod)      = (assemblyTuningPFcorr.(currPeriod).data - mean(assemblyTuningPFcorr.(currPeriod).ui, 2)) ./ std(assemblyTuningPFcorr.(currPeriod).ui, [], 2);
%     assemblyTuningSpatialInfo_zscore.(currPeriod) = (assemblyTuningSpatialInfo.(currPeriod).data - mean(assemblyTuningSpatialInfo.(currPeriod).ui, 2)) ./ std(assemblyTuningSpatialInfo.(currPeriod).ui, [], 2);
%     
    
end

% 
% fileName = [sessionName '.assemblyTunings_allPBEs_10ms.mat'];
% save(fullfile(storagePath, fileName), 'nPBEs', 'virtualOccupancy', ...
%     'assemblyTunings', 'assemblyTuningCorrMat', 'assemblyTuningPFcorr', 'assemblyTuningSpatialInfo') 

% 'assemblyTunings_zscore', 'assemblyTuningPFcorr_zscore', 'assemblyTuningSpatialInfo_zscore'


return

%% How assembly tunigns varies as a function of BD replay score

[assemblyTunings_sub, asTuningCorrMat_sub, asTuningPFcorr_sub, asTuningSpatialInfo_sub, ...
    assemblyTunings_sub_zscore, asTuningPFcorr_sub_zscore, asTuningSpatialInfo_sub_zscore] = deal(struct('pre', [], 'run', [], 'post', [])); % assembly tunings calculated based on PBEs with replay scores in different quartiles of the replay score distributions

nPBEs_subset = struct('pre', [], 'run', [], 'post', []);


% replayScoreMethods = {'rt_ts'; 'rt_ui'; 'rt_pf'; 'rt_ds'; ...
%                       'wc_ts'; 'wc_ui'; 'wc_pf'; 'wc_ds'};
% 
% replayScoreMethods_fullName = {'radon Integral - wPBE time swap'; ...
%                                'radon Integral - unit ID shuffle'; ...
%                                'radon Integral - place field shuffle'; ...
%                                'radon Integral - column cycle shuffle'; ...
%                                'weighted Corr - wPBE time swap'; ...
%                                'weighted Corr - unit ID shuffle'; ...
%                                'weighted Corr - place field shuffle'; ...
%                                'weighted Corr - column cycle shuffle'};
                  
replayScoreMethods = {'rt_ui'; 'wc_ui'; 'wc_ts'; 'wc_ds'; 'rt_ds'};
                      

replayScoreMethods_fullName = {'radon Integral - unit ID shuffle'; ...
                               'weighted Corr - unit ID shuffle' ; ...
                               'weighted Corr - wPBE time swap';
                               'weighted Corr - column cycle shuffle';
                               'radon Integral - column cycle shuffle'};
                           
fileName = [sessionName '.assemblyTuning_vs_replayScores_10ms.mat'];
                          
                           
for iperiod = 1:3
    
    
    currPeriod = periodNames{iperiod};
    currPBEs = PBEs.(currPeriod);

    posteriorProbMatrix = cell(nPBEs.(currPeriod), 1);
    for pbe = 1:nPBEs.(currPeriod) 
        posteriorProbMatrix{pbe} = currPBEs(pbe).posteriorProbMat;
    end
    posteriorProbMatrix = cell2mat(posteriorProbMatrix');


    virtualOccupancy.(currPeriod) = mean(posteriorProbMatrix, 2);
    virtualOccupancy.(currPeriod) = virtualOccupancy.(currPeriod) / sum(virtualOccupancy.(currPeriod)); % normalize 
    
end
                               
                           
nUIshuffles = 100;
nSeg        = 4;

nUIshuffles_perSeg = nUIshuffles/nSeg;


for irsm = 3%1:5%8
    
    replayScoreMethod = replayScoreMethods{irsm};
    replayScoreMethod_name = replayScoreMethods_fullName{irsm};
    
    fprintf(['/nReplay score method: ' replayScoreMethod_name ' ..'])
    
    
    for iperiod = 1:3 % behavioral epochs: pre, run, post
        
        currPeriod = periodNames{iperiod};
        
        fprintf(['/nepoch ' currPeriod ' ..'])
        
        
        
        
        currBDscores = [PBEs.(currPeriod).(replayScoreMethod)];

        nPBEs_subset.(currPeriod).(replayScoreMethod) = zeros(4,1);


        assemblyTunings_sub.(currPeriod).data.(replayScoreMethod)      = zeros(nUnits, nPosBins, 4);        
        assemblyTunings_sub.(currPeriod).ui.(replayScoreMethod)        = cell(4,1);

        asTuningSpatialInfo_sub.(currPeriod).data.(replayScoreMethod)  = zeros(nUnits, 4);
        asTuningPFcorr_sub.(currPeriod).data.(replayScoreMethod)       = zeros(nUnits, 4);
        asTuningCorrMat_sub.(currPeriod).data.(replayScoreMethod)      = zeros(nUnits, nUnits, 4);     


        asTuningSpatialInfo_sub.(currPeriod).ui.(replayScoreMethod)    = zeros(nUnits, 4, nUIshuffles);
        asTuningPFcorr_sub.(currPeriod).ui.(replayScoreMethod)         = zeros(nUnits, 4, nUIshuffles);


        for iPBEset = 1:4


            fprintf('\nCalculating assembly tunings based on PBEs in quartile %d of replay score ..', iPBEset)


            selectPBEs = find(currBDscores' > (iPBEset-1)*25 & currBDscores' <= iPBEset*25); 

            nPBEs_subset.(currPeriod).(replayScoreMethod)(iPBEset) = numel(selectPBEs);
            currbinnedPBEs = PBEs.(currPeriod);
            
            
            % calculate the assembly tunings based on the current PBEs
            
            ifShuffle = 0;
            currAssemblyTunings = calculateAssemblyTuningV7(currbinnedPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);

            assemblyTunings_sub.(currPeriod).data.(replayScoreMethod)(:, :, iPBEset) = currAssemblyTunings;


            meanFR  = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* currAssemblyTunings, 2) ./ sum(virtualOccupancy.(currPeriod));
            asTuningSpatialInfo_sub.(currPeriod).data.(replayScoreMethod)(:, iPBEset) = sum(repmat(virtualOccupancy.(currPeriod)', [nUnits 1]).* (currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(currAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);

            asTuningCorrMat_sub.(currPeriod).data.(replayScoreMethod)(:, :, iPBEset) = corr(permute(currAssemblyTunings, [2 1]), spatialTunings_merge');
            asTuningPFcorr_sub.(currPeriod).data.(replayScoreMethod)(:, iPBEset)     = diag(asTuningCorrMat_sub.(currPeriod).data.(replayScoreMethod)(:, :, iPBEset));

            
            

            % unit ID shuffle PBEs
            
            ifShuffle = 1;
            currAssemblyTunings = zeros(nUnits, nPosBins, nUIshuffles); % we are not going to store this variable, only need it for z-scoring 
            currPFcorr          = zeros(nUnits, nUIshuffles);
            currSpatialInfo     = zeros(nUnits, nUIshuffles);

            currOccupancy       = virtualOccupancy.(currPeriod);


            fprintf('\nCalculating assembly tunings based on unit ID shuffle surrogate PBEs ..')

            for seg = 1:nSeg
                seg
                parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)

                    shuffleAssemblyTunings = calculateAssemblyTuningV7(currbinnedPBEs, spikes, selectPBEs, clusterQuality, ifShuffle);
                    currAssemblyTunings(:, :, i_ui) = shuffleAssemblyTunings;

                    currCorrMat         = corr(shuffleAssemblyTunings', spatialTunings_merge'); 
                    currPFcorr(:, i_ui) = diag(currCorrMat);

                    meanFR                   = sum(repmat(currOccupancy', [nUnits 1]).* shuffleAssemblyTunings, 2) ./ sum(currOccupancy); % The calculation of mean firing rate could be a bit questionable
                    currSpatialInfo(:, i_ui) = sum(repmat(currOccupancy', [nUnits 1]).* (shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])) .* log2(shuffleAssemblyTunings ./ repmat(meanFR, [1 nPosBins])), 2);

                end
            end

            assemblyTunings_sub.(currPeriod).ui.(replayScoreMethod){iPBEset}           = currAssemblyTunings;
            asTuningPFcorr_sub.(currPeriod).ui.(replayScoreMethod)(:, iPBEset, :)      = permute(currPFcorr, [1 3 2]);
            asTuningSpatialInfo_sub.(currPeriod).ui.(replayScoreMethod)(:, iPBEset, :) = permute(currSpatialInfo, [1 3 2]);



            % comparing actual assembly tunings with the assembly tunings based on
            % unit ID shuffles
            assemblyTunings_sub_zscore.(currPeriod).(replayScoreMethod)(:, :, iPBEset) = (assemblyTunings_sub.(currPeriod).data.(replayScoreMethod)(:, :, iPBEset) - mean(currAssemblyTunings, 3)) ./ std(currAssemblyTunings, [], 3);

            asTuningPFcorr_sub_zscore.(currPeriod).(replayScoreMethod)(:, iPBEset)      = (asTuningPFcorr_sub.(currPeriod).data.(replayScoreMethod)(:, iPBEset) - mean(asTuningPFcorr_sub.(currPeriod).ui.(replayScoreMethod)(:, iPBEset, :), 3)) ./ std(asTuningPFcorr_sub.(currPeriod).ui.(replayScoreMethod)(:, iPBEset, :), [], 3);
            asTuningSpatialInfo_sub_zscore.(currPeriod).(replayScoreMethod)(:, iPBEset) = (asTuningSpatialInfo_sub.(currPeriod).data.(replayScoreMethod)(:, iPBEset) - mean(asTuningSpatialInfo_sub.(currPeriod).ui.(replayScoreMethod)(:, iPBEset, :), 3)) ./ std(asTuningSpatialInfo_sub.(currPeriod).ui.(replayScoreMethod)(:, iPBEset, :), [], 3);


        end

    end
    
    save(fullfile(storagePath, fileName), 'nPBEs_subset', ...
        'assemblyTunings_sub', 'asTuningCorrMat_sub', 'asTuningPFcorr_sub', 'asTuningSpatialInfo_sub', ...
        'assemblyTunings_sub_zscore', 'asTuningPFcorr_sub_zscore', 'asTuningSpatialInfo_sub_zscore')
end





end
