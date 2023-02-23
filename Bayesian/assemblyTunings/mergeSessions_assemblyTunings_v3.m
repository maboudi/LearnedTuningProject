
clear; clc

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

sub = dir(parentDir);
nSessions    = length(parentDir) - 2; % subtracting two becuase of '.' and '..'
sessionNames = cell(nSessions, 1);


epochNames = {'pre'; 'run'; 'post'};


replayScoreMethods = {'rt_ui'; 'wc_ui'; 'wc_ts'};


replayScoreMethods_fullName = {'radon Integral - unit ID shuffle'; ...
                               'weighted Corr - unit ID shuffle' ; ...
                               'weighted Corr - wPBE time swap'};


                           
% initialize 
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};


    % PF correlation
    assemblyTuning_pooled.(currEpoch).data = cell(nSessions, 1);
    assemblyTuning_pooled.(currEpoch).ui   = cell(nSessions, 1);
    
    
    assemblyTuningPFcorr_pooled.(currEpoch).data   = cell(nSessions, 1);
    assemblyTuningPFcorr_pooled.(currEpoch).ui     = cell(nSessions, 1);

    assemblyTuningPFcorr_zscore_pooled.(currEpoch) = cell(nSessions, 1);



    % PF correlation and spatial information across quartiles

    for irsm = 1:3 % index of replay score method

        replayScoreMethod = replayScoreMethods{irsm};

        % PF correlation
        asTuningPFcorr_sub_pooled.(currEpoch).data.(replayScoreMethod)   = cell(nSessions, 1); 
        asTuningPFcorr_sub_pooled.(currEpoch).ui.(replayScoreMethod)     = cell(nSessions, 1);

        asTuningPFcorr_sub_zscore_pooled.(currEpoch).(replayScoreMethod) = cell(nSessions, 1);

    end


end



for iSession = 3
    
    
    sessionName = sub(iSession+2).name
    
    
    fileName = [sessionName '.assemblyTunings_allPBEs.mat'];
    load(fullfile(parentDir, sessionName, 'assemblyTunings', fileName), 'nPBEs', 'virtualOccupancy', ...
            'assemblyTunings', 'assemblyTuningCorrMat', 'assemblyTuningPFcorr', 'assemblyTuningSpatialInfo', ...
            'assemblyTunings_zscore', 'assemblyTuningPFcorr_zscore', 'assemblyTuningSpatialInfo_zscore')

    fileName = [sessionName '.assemblyTuning_vs_replayScores.mat'];
    load(fullfile(parentDir, sessionName, 'assemblyTunings', fileName), 'nPBEs_subset', ...
            'assemblyTunings_sub', 'asTuningCorrMat_sub', 'asTuningPFcorr_sub', 'asTuningSpatialInfo_sub', ...
            'assemblyTunings_sub_zscore', 'asTuningPFcorr_sub_zscore', 'asTuningSpatialInfo_sub_zscore')
    
    
        
        
    for iEpoch = 1:3
    
        currEpoch = epochNames{iEpoch};
        
        
        % assembly Tunings
        assemblyTuning_pooled.(currEpoch).data{isession} = assemblyTunings.(currEpoch).data;
        assemblyTuning_pooled.(currEpoch).ui{isession}   = assemblyTunings.(currEpoch).ui;
        
        
        
        % PF correlation
        assemblyTuningPFcorr_pooled.(currEpoch).data{iSession}   = assemblyTuningPFcorr.(currEpoch).data;
        assemblyTuningPFcorr_pooled.(currEpoch).ui{iSession}     = assemblyTuningPFcorr.(currEpoch).ui;

        assemblyTuningPFcorr_zscore_pooled.(currEpoch){iSession} = assemblyTuningPFcorr_zscore.(currEpoch);


        % PF correlation and spatial information across quartiles

        for irsm = 1:3

            replayScoreMethod = replayScoreMethods{irsm};


            % PF correlation
            asTuningPFcorr_sub_pooled.(currEpoch).data.(replayScoreMethod){iSession}   = asTuningPFcorr_sub.(currEpoch).data.(replayScoreMethod); 
            asTuningPFcorr_sub_pooled.(currEpoch).ui.(replayScoreMethod){iSession}     = asTuningPFcorr_sub.(currEpoch).ui.(replayScoreMethod);

            asTuningPFcorr_sub_zscore_pooled.(currEpoch).(replayScoreMethod){iSession} = asTuningPFcorr_sub_zscore.(currEpoch).(replayScoreMethod);


        end
    
    end
    
end




%%

figure;

allCorrValues = [cell2mat(assemblyTuningPFcorr_zscore_pooled.pre); cell2mat(assemblyTuningPFcorr_zscore_pooled.run); cell2mat(assemblyTuningPFcorr_zscore_pooled.post)];
bins = linspace(min(allCorrValues), max(allCorrValues)+0.01, 15);



colors = {'b'; 'none';'r'};
edgeColors = {'none'; 'k'; 'none'};
for iEpoch = 1:3
    
    currPeriod = epochNames{iEpoch};
    
    h = histc(cell2mat(assemblyTuningPFcorr_zscore_pooled.(currPeriod)), bins); h(end) = [];
    a(iEpoch ) = bar(bins(1:end-1)+diff(bins(1:2))/2, h/sum(h), 'Edgecolor', edgeColors{iEpoch}, 'facecolor', colors{iEpoch}, 'facealpha', 0.3);
    hold on
    
    
    rr = ylim;
    medianCorr = nanmedian(cell2mat(assemblyTuningPFcorr_zscore_pooled.(currPeriod)));
    
    
%     shuffle = cell2mat(assemblyTuningPFcorr_pooled.(currPeriod).ui);
%     medianCorr_shuffle = nanmedian(shuffle(:));
    
    
    cl = colors{iEpoch}; 
    if iEpoch == 2; cl = 'k'; end
 
%     line([medianCorr_shuffle medianCorr_shuffle], [0 rr(2)], 'color', cl, 'linewidth', 2, 'linestyle', ':')
    line([medianCorr medianCorr], [0 rr(2)], 'color', cl, 'linewidth', 2)
    
end

% xlim([-1 1])

% set(gca,'YColor','none')
% xticklabels([])

rr = legend(a, 'PRE', 'RUN', 'POST', 'box', 'off');
rr.FontSize = 14;

set(gca, 'box', 'off', 'fontsize', 12, 'linewidth', 2)

xlabel({'correlation bw assembly tunings and place fields'}, 'fontsize', 14)
ylabel('fraction of units', 'fontsize', 14)
title('all datasets', 'fontsize', 14)

set(gca, 'box', 'off')



%%   Assembly tunings versus BD replay scores     
        

replayScoreMethodNumber = 1;


% correlation of assembly tunings with place fields

figure;
set(gcf, 'position', [3000 50 410 900])
colors = {'b'; 'k'; 'r'};
edgeColors = {'k', 'g', 'k'};
panelNumbers = [3 4 5; 7 8 9; 11 12 13];


replayScoreMethod = replayScoreMethods{replayScoreMethodNumber};
replayScoreMethod_fullName = replayScoreMethods_fullName{replayScoreMethodNumber};



for iEpoch = 1:3

    currPeriod = epochNames{iEpoch};

    subplot(14,1, panelNumbers(iEpoch, :))
    set(gca, 'fontsize', 12, 'linewidth', 2)
    
    h = violinplot(cell2mat(asTuningPFcorr_sub_pooled.(currPeriod).data.(replayScoreMethod)));

    for ii = 1:4
        h(ii).ViolinColor = colors{iEpoch};
        h(ii).EdgeColor   = edgeColors{iEpoch};
        h(ii).ShowData    = 0;
    end

    if iEpoch == 3
        xticklabels({'0-25%'; '25-50%'; '50-75%'; '75-100%'})
    else
        xticklabels([])
    end

    
    if iEpoch == 2
        ylabel('correlation bw assembly tunings and place fields', 'fontsize', 18)
    end
    if iEpoch == 3; xlabel('Bayesian replay scores', 'fontsize', 14); end
    xtickangle(45)


    
    ylim([-1 1])
    
    grid on


    
    
    hold on
    
    data = cell2mat(asTuningPFcorr_sub_pooled.(currPeriod).data.(replayScoreMethod));
    shuffle = cell2mat(asTuningPFcorr_sub_pooled.(currPeriod).ui.(replayScoreMethod));
    
    for ii = 1:4
        
        temp = shuffle(:, ii, :);
        corrWithPFs_ui = temp(:);
        
        dataMedian    = nanmedian(data(:, ii));
        shuffleMedian = nanmedian(corrWithPFs_ui);
        line([ii-0.3 ii+0.3], [shuffleMedian shuffleMedian], 'color', colors{iEpoch}, 'linestyle', ':', 'linewidth', 2)

        pval = ranksum(data(:, ii), corrWithPFs_ui, 'tail', 'right');
        
        if pval < 0.001
            text(ii-0.5, -0.5, sprintf('p<0.001'), 'color', [0.1 0.1 0.1], 'fontsize', 8)
        else
            text(ii-0.5, -0.5, sprintf('p = %.3f', pval), 'color', [0.1 0.1 0.1], 'fontsize', 8)
        end
% %         text(ii-0.5, dataMedian+0.4, sprintf('n = %d', nPBEs_subset.(currPeriod)(ii)), 'color', 'b', 'fontsize', 8)
    end

    hold off

    tt = ylim;
    text(0.5, tt(2), currPeriod, 'fontsize', 12, 'fontweight', 'bold', 'horizontalAlignment', 'center')

    dim = [0.1 0.8 0.5 .1];
    annotation('textbox',dim,'String', {'all datasets'; ['replay scoring method: ' replayScoreMethod_fullName]}, ...
          'FitBoxToText','on', 'Interpreter', 'none'); % sessionName;

end


return



%% plot assembly tunings for an individual session
% (THIS PART OF CODE SHOULD BE IN ANOTHER M.FILE)

sessionNumber = 4;

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/concat_GreatLakes_datasets_temp/';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name;

basePath = fullfile(parentDir, sessionName);

% load(fullfile(basePath, 'PopulationBurstEvents', [sessionName '.PBEInfo.mat']), 'PBEInfo', 'acceptedIdx')
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr')
load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']))
load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores.mat']))

spikes = spikes_pyr;

% for better plots sort the units based on their place field peak position
nUnits = numel(spikes);
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
spatialTunings_merge = zeros(nUnits, nPosBins);

peakPosBins = zeros(nUnits, 1);
for iUnit = 1:nUnits
    peakPosBins(iUnit) = spikes(iUnit).peakPosBin.uni;
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
end
    
[~, pf_sortIndices] = sort(peakPosBins, 'descend');





figure; 

epochNames = {'pre'; 'run'; 'post'};

for ii = 1:3
    
    currPeriod = epochNames{ii};
    ax(ii) = subplot(1,3,ii);
    
    
    currAT = assemblyTunings_sub.(currPeriod).data.wc_ts(pf_sortIndices, :, 4);
%     currAT = assemblyTunings.(currPeriod).data(pf_sortIndices, :);
    currAT = currAT ./ repmat(max(currAT, [], 2), [1 size(currAT, 2)]);
    
    [nUnits, nPosBins] = size(currAT);
    
    imagesc(2*(1:nPosBins), 1:nUnits, currAT)
    colormap(ax(ii), 'jet')
    xlabel('track position', 'fontsize', 12); 
    if ii == 1 
        ylabel('unit', 'fontsize', 12)
    end
    title([currPeriod sprintf('(n=%d)', nPBEs.(currPeriod))]);
%     caxis([0 1])
end



figure; 

spatialTunings = spatialTunings_merge(pf_sortIndices, :);
spatialTunings = spatialTunings ./ repmat(max(spatialTunings, [], 2), [1 size(currAT, 2)]);

    
imagesc(2*(1:nPosBins), 1:nUnits, spatialTunings)
colormap('jet')
xlabel('track position', 'fontsize', 12); 
if ii == 1 
    ylabel('unit', 'fontsize', 12)
end
title([currPeriod sprintf('(n=%d)', nPBEs.(currPeriod))]);
caxis([0 1])



