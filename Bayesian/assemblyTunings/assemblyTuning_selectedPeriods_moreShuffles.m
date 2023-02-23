function assemblyTuning_selectedPeriods_moreShuffles(sessionNumber)

% we might need to consider only stable units when calculating the assembly
% tuning for a specific unit


sz = getenv('SLURM_CPUS_PER_TASK');

theCluster = parcluster('local');

JobFolder = sprintf('/home/kmaboudi/.matlab/trash/job_%s', sessionNumber);
mkdir(JobFolder)
theCluster.JobStorageLocation = JobFolder;

p = parpool(theCluster, str2double(sz));


addpath(genpath('/nfs/turbo/umms-kdiba/NCMLproject/ReplayPreplayAnalyses'))

parentDir = '/nfs/turbo/umms-kdiba/Kourosh/NCMLproject/concat_GreatLakes_datasets_temp';
% parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name


basePath    = fullfile(parentDir, sessionName);
storagePath = fullfile(basePath, 'assemblyTunings');
mkdir(storagePath)

load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')


% loading the brain state info if it exists (I calculated the brains states for Grosmark ... 
% and Bapun's sessions from scratch but for Hiro's data I trusted whatever he calculated)

fileName = fullfile(basePath, 'spikes', [sessionName '.brainState.mat']);
if isfile(fileName)
    load(fileName, 'brainState')
else 
    brainState = fileInfo.brainStates;
end

load(fullfile(basePath, 'PopulationBurstEvents', [sessionName '.PBEInfo.mat']), 'sdat')



% behavioral epochs

behavior = fileInfo.behavior;

startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
startT.post = behavior.time(2,2); endT.post = startT.post + 5*60*60; %behavior.time(3,2); 


% clusterQuality = clusterQuality; 

spikes   = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
nUnits   = numel(spikes); 



% non-directional spatial tuning

spatialTunings = nan(nUnits, nPosBins);
peakPFfiring   = nan(nUnits, 1);
PFpeakLoc        = nan(nUnits, 1);

for iUnit = 1:nUnits
    spatialTunings(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni; 
    [~, PFpeakLoc(iUnit)] = max(spatialTunings(iUnit, :));
    
    peakPFfiring(iUnit) = spikes(iUnit).peakFR.uni;

end
[~, PFSortIdx] = sort(PFpeakLoc, 'ascend');


% NREM and QW periods with MUA(z) below zero, also REM periods

sdat_t = (1:length(sdat))*1e-3;

minLowFRDur = 0.1;
lowFRthresh = 0;
minREMDur    = 0.5;
min_interREM = 0.5;

if ismember(sessionNumber, [1:5 7:9])

    theta_intrp = interp1(brainState.thetaTimePnts, brainState.theratio, sdat_t);
    sw_intrp    = interp1(brainState.sw_emg_timePnts, brainState.slowWave, sdat_t);
    emg_intrp   = interp1(brainState.sw_emg_timePnts, brainState.emg, sdat_t);

    
    % NREM non-PBE periods
    lowFRIdx = sdat < lowFRthresh & (sw_intrp' > brainState.swthresh | theta_intrp' < 1);
    
    % REM periods
    REMIdx = sw_intrp' < brainState.swthresh & theta_intrp' > 1 & emg_intrp' < brainState.emgThresh;

    
elseif ismember(sessionNumber, [6 10:15])
    
    brainState(:,2) = brainState(:, 2) - 1e-3;
    brainState(:,4) = brainState(:, 3);
    
    timePnts = brainState(:, 1:2)';
    timePnts = timePnts(:);
    
    stateIdx = brainState(:, 3:4)';
    stateIdx = stateIdx(:);
    
    stateIdx_intrp = interp1(timePnts, stateIdx, sdat_t);
        
    % NREM non-PBE periods    
    lowFRIdx = sdat < lowFRthresh & (ismember(stateIdx_intrp, [1 3]))';

    % REM peridos
    REMIdx = stateIdx_intrp == 2;
    REMIdx = REMIdx';
    
end



% NREM bouts
try 
    crossup   = find([0; diff(lowFRIdx)] == 1);
    crossdown = find([0; diff(lowFRIdx)] == -1);
catch
    crossup   = find([0 diff(lowFRIdx)] == 1);
    crossdown = find([0 diff(lowFRIdx)] == -1);
end


if crossdown(1) < crossup(1); crossdown(1) = []; end
if crossup(end) > crossdown(end); crossup(end) = []; end

lowFRBouts = sdat_t([crossup crossdown]); 
lowFRBouts(diff(lowFRBouts')' < minLowFRDur, :) = [];



% REM bouts
try 
    crossup   = find([0; diff(REMIdx)] == 1);
    crossdown = find([0; diff(REMIdx)] == -1);
catch 
    crossup   = find([0 diff(REMIdx)] == 1);
    crossdown = find([0 diff(REMIdx)] == -1);
end


if crossdown(1) < crossup(1); crossdown(1) = []; end
if crossup(end) > crossdown(end); crossup(end) = []; end

REMBouts = sdat_t([crossup crossdown]);


% concatenate REM bouts that are less than min_interREM apart
% We don't do it for NREM periods as the concatenation might lead to
% inclusion of population burst events .. 
    
tmp_rem = REMBouts(1,:);
finalREMBouts = [];

for ii = 2:size(REMBouts,1)
    if (REMBouts(ii,1)-tmp_rem(2)) < min_interREM
        tmp_rem = [tmp_rem(1) REMBouts(ii, 2)]; % merge two PBEs
    else
        finalREMBouts = [finalREMBouts; tmp_rem];
        tmp_rem = REMBouts(ii,:);
    end
end

finalREMBouts = [finalREMBouts;tmp_rem]; 

REMBouts = finalREMBouts;

REMdur  = REMBouts(:, 2) - REMBouts(:, 1);
% REMdist = REMBouts(2:end, 1) - REMBouts(1:end-1, 2); 

REMBouts(REMdur < minREMDur, :) = [];



%% assembly Tunings 
% Spatial information and correlation with place fields of assesmbly tunings
% and corresponding z-scores by comparing the values againts distibutions calculated based on unit ID shuffle surrogates


bouts = lowFRBouts; 

% binDur = 0.02;

epochNames = {'pre'; 'run'; 'post'};
dataTypes  = {'data'; 'ui'};

nUIshuffles = 100;
nSeg        = 1;

nUIshuffles_perSeg = nUIshuffles/nSeg;

for binDur = 0.02

    
    for iEpoch = [1 3]

        currEpoch = epochNames{iEpoch};
        epoch     = [startT.(currEpoch) endT.(currEpoch)];
        

        currEpoch_includeIdx = bouts(:, 1) > epoch(1) & bouts(:, 1) < epoch(2);
        currBouts = bouts(currEpoch_includeIdx, :);


        fprintf(['\nProcessing ' currEpoch ' ,binDur = ' num2str(binDur) ' ..'])
    
        % actual data

        fprintf('\nCalculating the assembly tunings based on actual PBEs ..')
        

        ifShuffle = 0;
        assemblyTunings.(currEpoch).data = calculateAssemblyTuning_selectedPeriods(currBouts, spikes, binDur, clusterQuality, ifShuffle);    

%         [assemblyTunings.(currEpoch).data, concatBinCenters.(currEpoch), ...
%                 concatBoutIDs.(currEpoch), boutFirings_unit.(currEpoch).data] = calculateAssemblyTuning_selectedPeriods(currBouts, spikes, binDur, clusterQuality, ifShuffle);

        
        fprintf('\nCalculating the assembly tunings based on unit identity shuffled PBEs ..')
        ifShuffle = 1;
        
        

        currAssemblyTunings = nan(nUnits, nPosBins, nUIshuffles);

        for seg = 1:nSeg

            parfor i_ui = ((seg-1)*nUIshuffles_perSeg +1):(seg*nUIshuffles_perSeg)
                currAssemblyTunings(:, :, i_ui) = calculateAssemblyTuning_selectedPeriods(currBouts, spikes, binDur, clusterQuality, ifShuffle);
            end
        end

        assemblyTunings.(currEpoch).ui = currAssemblyTunings;



        for idata = 1:2
            
            currDataType = dataTypes{idata};
            
            currAssemblyTunings = assemblyTunings.(currEpoch).(currDataType);
            nSamples = size(currAssemblyTunings, 3);

            assemblyTuningPFcorr.(currEpoch).(currDataType)         = nan(nUnits, nSamples);
            assemblyTuningPFKLdiv.(currEpoch).(currDataType)        = nan(nUnits, nSamples);
            assemblyTuningPF_prefPosDist.(currEpoch).(currDataType) = nan(nUnits, nSamples);
            assemblyTuning_prefPos.(currEpoch).(currDataType)       = nan(nUnits, nSamples);

            for is = 1:nSamples

                allCorrs = corr(currAssemblyTunings(:, :, is)', spatialTunings');
                assemblyTuningPFcorr.(currEpoch).(currDataType)(:, is) = diag(allCorrs);
    
                allKLdiv = calKLDivergence(currAssemblyTunings(:, :, is)', spatialTunings');
                assemblyTuningPFKLdiv.(currEpoch).(currDataType)(:, is) = diag(allKLdiv);
                
                [assemblyTuningPF_prefPosDist.(currEpoch).(currDataType)(:, is), assemblyTuning_prefPos.(currEpoch).(currDataType)(:, is)] = ...
                    calDistPrefPosition(currAssemblyTunings(:, :, is), PFpeakLoc, []);
            end

        end
        
        assemblyTuningPFcorr_pscore.(currEpoch)         = calPrctile(assemblyTuningPFcorr.(currEpoch));
        assemblyTuningPFKLdiv_pscore.(currEpoch)        = calPrctile(assemblyTuningPFKLdiv.(currEpoch));
        assemblyTuningPF_prefPosDist_pscore.(currEpoch) = calPrctile(assemblyTuningPF_prefPosDist.(currEpoch));
    end


    fileName = sprintf('%s.assemblyTunings_NREM_nonPBE_%.3f.mat', sessionName, binDur);

    save(fullfile(storagePath, fileName), ...
        'lowFRthresh', ...
        'lowFRBouts', ...
        'assemblyTunings', ...
        'assemblyTuning_prefPos', ...
        'assemblyTuningPF_prefPosDist', ...
        'assemblyTuningPF_prefPosDist_pscore', ...
        'assemblyTuningPFcorr', ...
        'assemblyTuningPFcorr_pscore', ...
        'assemblyTuningPFKLdiv', ...
        'assemblyTuningPFKLdiv_pscore');

    
end


return


%% 

% binDur  = 4*60*60; % 15 minutes
% stepDur = 600; % 5 minutes

for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
    
    fprintf(['\nProcessing ' currEpoch ' ..'])
    
%     totalT = endT.(currEpoch) - startT.(currEpoch);
%     nBins  = floor((totalT - binDur)/stepDur) + 1;
     
%     binStarts   = startT.(currEpoch) + (0:nBins-1)*stepDur;
%     binEnds     = binStarts + binDur;
    
    nBins       = 1;
    binStarts   = startT.(currEpoch);
    binEnds     = endT.(currEpoch);

%     binCenters.(currEpoch)  = binStarts + binDur/2;
    
    
    for idata = 1%:2
        
        currDataType = dataTypes{idata};
        
        learnedTunings2 = nan(nUnits, nPosBins, nBins);   
        
        unit_fr.(currEpoch).(currDataType)               = nan(nUnits, nBins);
        spatialBinCorrelation.(currEpoch).(currDataType) = nan(nUnits, nBins);
%         spatialBinCorrPval.(currEpoch).(currDataType)    = nan(nUnits, nBins);
        KLdivergence.(currEpoch).(currDataType)          = nan(nUnits, nBins);
        KSstatistic.(currEpoch).(currDataType)           = nan(nUnits, nBins);

        for iBin = 1%:nBins

            idx = concatBinCenters.(currEpoch) > binStarts(iBin) & concatBinCenters.(currEpoch) < binEnds(iBin);
            
            if ~all(~idx)
                currLearnedTunings = assemblyTunings.(currEpoch).(currDataType)(:, :, idx);
                currLearnedTunings = nanmean(currLearnedTunings .* repmat(permute(boutFirings_unit.(currEpoch).(currDataType)(:, idx), [1 3 2]), [1 nPosBins]), 3);
                learnedTunings2(:, :, iBin) = currLearnedTunings ./ repmat(sum(currLearnedTunings, 2), [1 nPosBins]);
                learnedTunings2(:, :, iBin) = learnedTunings2(:, :, iBin) + eps;


                for iUnit = 1:nUnits

                    unit_fr.(currEpoch).(currDataType)(iUnit, iBin) = nanmean(boutFirings_unit.(currEpoch).data(iUnit, idx)) * 0.05;

                    spatialBinCorrelation.(currEpoch).(currDataType)(iUnit, iBin) = corr(learnedTunings2(iUnit, :, iBin)', spatialTunings(iUnit, :)', 'tail', 'right');

                    KLdivergence.(currEpoch).(currDataType)(iUnit, iBin) = sum(learnedTunings2(iUnit, :, iBin) .* (log(learnedTunings2(iUnit, :, iBin)) - log(spatialTunings(iUnit, :))));  
                    KSstatistic.(currEpoch).(currDataType)(iUnit, iBin)  = max(abs(cumsum(learnedTunings2(iUnit, :, iBin)) - cumsum(spatialTunings(iUnit, :))));

                end
            end 
        end
        
        learnedTunings.(currEpoch).(currDataType) = learnedTunings2; 
    end
end


%% plots

placeCells = find(peakPFfiring >= 2);

for iUnit = 1:nUnits

    
if ~ismember(iUnit, placeCells)
    continue
end

plotheight = 10;
plotwidth  = 15;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



ax(1) = subplot(14,1,2);
hold on
for iEpoch = [1 3]
    currEpoch = epochNames{iEpoch};

    plot(binCenters.(currEpoch), unit_fr.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7]);
end
xlim([0 binCenters.post(end)])
xticklabels([])


set(gca, 'fontsize', fontsize)


yl = ylim(ax(1));

text(startT.run, yl(1)+2*diff(yl), 'run', 'fontsize', 10)
text(mean([startT.pre endT.pre]), yl(1)+2*diff(yl), 'pre', 'fontsize', 10)
text(mean([startT.post endT.post]), yl(1)+2*diff(yl), 'post', 'fontsize', 10)

text(startT.pre, yl(1)+3*diff(yl),  sprintf('unit-%d', iUnit), 'fontweight', 'bold')


ax(2) = subplot(14,1, 3:7);
hold on

currSpatialTuning = spikes(iUnit).spatialTuning_smoothed.uni;

for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
    
    learnedTuning2 = squeeze(learnedTunings.(currEpoch).data(iUnit, :, :, 1));
    learnedTuning2 = learnedTuning2 ./ repmat(max(learnedTuning2, [], 1), [size(learnedTuning2, 1) 1]);
    
    imagesc(binCenters.(currEpoch), 1:nPosBins, learnedTuning2); colormap('jet')

end

plot(ax(2), 500*currSpatialTuning + startT.run - diff([startT.run endT.run])/2 , 1:numel(currSpatialTuning), 'linewidth', 2, 'color', [0 0 0 0.5], 'DisplayName', 'track spatial tuning')


yl = ylim(ax(2));

line(500*[0 1] + startT.run - diff([startT.run endT.run])/2 , [yl(2)+5 yl(2)+5], 'linewidth', 2, 'color', [0 0 0 0.5])
text(startT.run - diff([startT.run endT.run])/2, yl(2)+15, '1Hz', 'fontsize', 6, 'color', 'k')

ylim([0 nPosBins+10])


[maxFR, peakPos] = max(currSpatialTuning);

b1 = find(currSpatialTuning < maxFR * 0.25 & (1:nPosBins)' < peakPos, 1, 'last');
b2 = find(currSpatialTuning < maxFR * 0.25 & (1:nPosBins)' > peakPos, 1, 'first');



xlim([0 binCenters.post(end)])
xl = xlim;
try
line(xl, [b1, b1], 'color', 'w', 'linestyle', '--', 'linewidth', 0.5)
catch
end

try
line(xl, [b2, b2], 'color', 'w', 'linestyle', '--', 'linewidth', 0.5)
catch
end

xticklabels([])

set(gca, 'fontsize', fontsize)



% % spatial bin correlation between the learned tunings and track spatial
% bin tuning

ax(3) = subplot(14, 1, 8:9);
hold on

for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
    plot(binCenters.(currEpoch), spatialBinCorrelation.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'avg of multiple time bins');
%     plot(binCenters.(currEpoch), spatialBinCorrelation.(currEpoch).ui(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.3], 'DisplayName', 'avg of multiple time bins');

    
%     pvalues = spatialBinCorrPval.(currEpoch).data(iUnit, :);

    nBins = numel(binCenters.(currEpoch));
    pvalues = nan(nBins, 1);
    for iBin = 1:nBins
        try
            pvalues(iBin) = ranksum(spatialBinCorrelation.(currEpoch).data(iUnit, iBin), spatialBinCorrelation.(currEpoch).ui(:), 'tail', 'right');
        catch
        end
    end
    
    significance = 1 - pvalues;
%     significance = nan(size(pvalues));
%     significance(pvalues < 0.05) = 1;
    
%     plot(binCenters.(currEpoch), significance, '.', 'markersize', 5, 'color', [.7 .7 .7 .7]);
    plot(binCenters.(currEpoch), significance, 'linewidth', 1, 'color', [.7 .7 .7 .7]);

    
    xl = binCenters.(currEpoch)([1 end]);
    
    line(xl, [nanmedian(spatialBinCorrelation.(currEpoch).ui(:)) nanmedian(spatialBinCorrelation.(currEpoch).ui(:))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median_indiv bins');
    line(xl, [nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 25) nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
    line(xl, [nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 75) nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')


%     line(xl, [nanmedian(spatialBinCorrelation.(currEpoch).data(:)) nanmedian(spatialBinCorrelation.(currEpoch).data(:))], 'color', [0 0 0 0.5], 'linewidth', 1, 'DisplayName', 'pooled ui median_indiv bins');
%     line(xl, [nanprctile(spatialBinCorrelation.(currEpoch).data(:), 25) nanprctile(spatialBinCorrelation.(currEpoch).data(:), 25)], 'color', [0 0 0 0.5], 'linewidth', 0.5, 'linestyle', '--')
%     line(xl, [nanprctile(spatialBinCorrelation.(currEpoch).data(:), 75) nanprctile(spatialBinCorrelation.(currEpoch).data(:), 75)], 'color', [0 0 0 0.5], 'linewidth', 0.5, 'linestyle', '--')


    
end
xticklabels([])
xlim([0 binCenters.post(end)])
% xl = xlim;

% line(xl, [nanmedian(spatialBinCorrelation.(currEpoch).ui(:)) nanmedian(spatialBinCorrelation.(currEpoch).ui(:))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median_indiv bins');
% line(xl, [nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 25) nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
% line(xl, [nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 75) nanprctile(spatialBinCorrelation.(currEpoch).ui(:), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')

set(gca, 'fontsize', fontsize)



% % KL distance: bits of infomration lost when I use learned tunings to
% describ the information encoded in the track spatail tunings
% 
% ax(4) = subplot(14, 1, 10:11);
% hold on
% 
% for iEpoch = [1 3]
%     
%     currEpoch = epochNames{iEpoch};
% 
%     plot(binCenters.(currEpoch), KLdivergence.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'avg of multiple time bins');
% 
% end
% 
% xticklabels([])
% xlim([0 binCenters.post(end)])
% xl = xlim;
% 
% line(xl, [nanmedian(KLdivergence.(currEpoch).ui(iUnit, :)) nanmedian(KLdivergence.(currEpoch).ui(iUnit, :))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median: 20msbin-resolved');
% line(xl, [nanprctile(KLdivergence.(currEpoch).ui(iUnit, :), 25) nanprctile(KLdivergence.(currEpoch).ui(iUnit, :), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
% line(xl, [nanprctile(KLdivergence.(currEpoch).ui(iUnit, :), 75) nanprctile(KLdivergence.(currEpoch).ui(iUnit, :), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
% 
% set(gca, 'fontsize', fontsize)
% 
% 


% % KS stats: the max difference between two CDFs
% 
% ax(5) = subplot(14, 1, 12:13);
% hold on
% 
% for iEpoch = [1 3]
%     
%     currEpoch = epochNames{iEpoch};
% 
%     plot(binCenters.(currEpoch), KSstatistic.(currEpoch).data(iUnit, :), 'linewidth', 2, 'color', [0 0 0 0.7], 'DisplayName', 'avg of multiple time bins');
% 
% end
% 
% 
% xlim([0 binCenters.post(end)])
% xl = xlim;
% 
% line(xl, [nanmedian(KSstatistic.(currEpoch).ui(iUnit, :)) nanmedian(KSstatistic.(currEpoch).ui(iUnit, :))], 'color', [0 0 0 0.3], 'linewidth', 1, 'DisplayName', 'pooled ui median: 20msbin-resolved');
% line(xl, [nanprctile(KSstatistic.(currEpoch).ui(iUnit, :), 25) nanprctile(KSstatistic.(currEpoch).ui(iUnit, :), 25)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
% line(xl, [nanprctile(KSstatistic.(currEpoch).ui(iUnit, :), 75) nanprctile(KSstatistic.(currEpoch).ui(iUnit, :), 75)], 'color', [0 0 0 0.3], 'linewidth', 0.5, 'linestyle', '--')
% 
% set(gca, 'fontsize', fontsize)
% 






linkaxes(ax, 'x')

ylabel(ax(1), 'PBE fr(Hz)', 'fontsize', fontsize)
ylabel(ax(2), 'position bin', 'fontsize', fontsize)
ylabel(ax(3), 'correlation', 'fontsize', fontsize)

xlabel(ax(3), 'time(sec)', 'fontsize', fontsize)


if iUnit < 10
    fileName = sprintf('%s_unit%d%d.pdf', sessionName, 0, iUnit);
else
    fileName = sprintf('%s_unit%d.pdf', sessionName, iUnit);
end

print(gcf, fileName, '-dpdf', '-painters')

close all



end

end
%% subfunctions

function output = nanprctile(x, percent)

output = prctile(x(~isnan(x)), percent);

end


function KLdiv = calKLDivergence(learnedTunings, spatialTunings)

% the units should be in columns of learned tunings or spatialTunings
% matrix
% 
% if size(spatialTunings, 2) > size(spatialTunings, 1)
%     spatialTunings = spatialTunings';
% end

spatialTunings = spatialTunings ./ repmat(sum(spatialTunings, 1), [size(spatialTunings, 1) 1]);
spatialTunings = spatialTunings + eps;

[~, nUnits] = size(spatialTunings);

% if size(learnedTunings, 1) ~= nPosBins
%     learnedTunings = learnedTunings';
% end

learnedTunings = learnedTunings ./ repmat(sum(learnedTunings, 1), [size(learnedTunings, 1) 1]);
learnedTunings = learnedTunings + eps;


KLdiv = nan(nUnits);

for iUnit = 1:nUnits
    currSpatialTuning = spatialTunings(:, iUnit);
    
    spatialTuningTerm = repmat(currSpatialTuning, [1 nUnits]);
    KLdiv(:, iUnit) = sum(spatialTuningTerm .* (log(spatialTuningTerm) - log(learnedTunings)), 1);
end

end




function [prefPosDist, asTuning_prefPos] = calDistPrefPosition(assemblyTunings, PF_prefPos, activeUnits)  
    
    nUnits = numel(PF_prefPos);
    nSets = size(assemblyTunings, 3);
    nPosBins = size(assemblyTunings, 2);
    
    if isempty(activeUnits)
        activeUnits = 1:nUnits;
    end
    
    prefPosDist   = nan(nUnits, nSets);
    asTuning_prefPos = nan(nUnits, nSets);

    for n = 1:nSets
       currAsTunings = assemblyTunings(:,:, n);
       [~, asTuning_prefPos(:, n)] = max(currAsTunings, [], 2);
       prefPosDist(activeUnits, n) = abs(asTuning_prefPos(activeUnits, n) - PF_prefPos(activeUnits));
       prefPosDist(activeUnits, n) = prefPosDist(activeUnits, n)/nPosBins;
    end

end


function prctileScore = calPrctile(variable)


[nUnits, nShuffles] = size(variable.ui);
prctileScore = zeros(nUnits, 1);

for iUnit = 1:nUnits
    
    if isnan(variable.data(iUnit))
        prctileScore(iUnit) = nan;
    else
        prctileScore(iUnit) = numel(find(variable.ui(iUnit, :) <= variable.data(iUnit)))/nShuffles * 100;
    end

end

end