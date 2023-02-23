
clear; clc

currSession = '/home/kouroshmaboudi/Documents/HMM_project/assemblyTuning_Feb12th/';
cd(currSession)

currDir = dir(currSession);
nSessions    = length(currDir) - 2; % subtracting two becuase of '.' and '..'
sessionNames = cell(nSessions, 1);

periodNames = {'PRE'; 'RUN'; 'POST'};
shuffles    = {'wPBEtimeswap'; 'unitIDshuffle'};
seqMetrics  = {'replayScore'; 'weightedCorr'};



for iperiod = 1:3

    currentPeriod = periodNames{iperiod};


    % PF correlation
    assemblyTuningPFcorr_pooled.(currentPeriod).data   = cell(nSessions, 1);
    assemblyTuningPFcorr_pooled.(currentPeriod).ui     = cell(nSessions, 1);

    assemblyTuningPFcorr_zscore_pooled.(currentPeriod) = cell(nSessions, 1);



    % Spatial information
    assemblyTuningSpatialInfo_pooled.(currentPeriod).data = cell(nSessions, 1);
    assemblyTuningSpatialInfo_pooled.(currentPeriod).ui   = cell(nSessions, 1);

    assemblyTuningSpatialInfo_zscore_pooled.(currentPeriod) = cell(nSessions, 1);



    % PF correlation and spatial information across quartiles

    for ishuffle = 1:2
        for iseq = 1:2

            currShuffleMethod = shuffles{ishuffle};
            currSeqMetric     = seqMetrics{iseq};

            % PF correlation
            asTuningPFcorr_sub_pooled.(currentPeriod).data.(currShuffleMethod).(currSeqMetric)   = cell(nSessions, 1); 
            asTuningPFcorr_sub_pooled.(currentPeriod).ui.(currShuffleMethod).(currSeqMetric)     = cell(nSessions, 1);
            
            asTuningPFcorr_sub_zscore_pooled.(currentPeriod).(currShuffleMethod).(currSeqMetric) = cell(nSessions, 1);

            
            
            
            % Spatial information
            asTuningSpatialInfo_sub_pooled.(currentPeriod).data.(currShuffleMethod).(currSeqMetric) = cell(nSessions, 1);
            asTuningSpatialInfo_sub_pooled.(currentPeriod).ui.(currShuffleMethod).(currSeqMetric)   = cell(nSessions, 1);
            
            asTuningSpatialInfo_sub_zscore_pooled.(currentPeriod).(currShuffleMethod).(currSeqMetric) = cell(nSessions, 1);

        end
    end


end



for iSession = [1:5 6 9 10]
    
    
    
    sessionName = currDir(iSession+2).name
    sessionNames{iSession} = sessionName;
    
    load(fullfile(currSession, sessionName, 'assemblyTunings_allPBEs'), ...
                'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore', ...
                'assemblyTuningSpatialInfo', 'assemblyTuningSpatialInfo_zscore')

    load(fullfile(currSession, sessionName, 'assemblyTuning_vs_replayScores.mat'), ...
                'asTuningPFcorr_sub', 'asTuningPFcorr_sub_zscore', ...
                'asTuningSpatialInfo_sub', 'asTuningSpatialInfo_sub_zscore')
    
    for iperiod = 1:3
    
        currentPeriod = periodNames{iperiod};
        % PF correlation
        assemblyTuningPFcorr_pooled.(currentPeriod).data{iSession}   = assemblyTuningPFcorr.(currentPeriod).data;
        assemblyTuningPFcorr_pooled.(currentPeriod).ui{iSession}     = assemblyTuningPFcorr.(currentPeriod).ui;

        assemblyTuningPFcorr_zscore_pooled.(currentPeriod){iSession} = assemblyTuningPFcorr_zscore.(currentPeriod);


        % Spatial information
        assemblyTuningSpatialInfo_pooled.(currentPeriod).data{iSession} = assemblyTuningSpatialInfo.(currentPeriod).data;
        assemblyTuningSpatialInfo_pooled.(currentPeriod).ui{iSession}   = assemblyTuningSpatialInfo.(currentPeriod).ui;

        assemblyTuningSpatialInfo_zscore_pooled.(currentPeriod){iSession} = assemblyTuningSpatialInfo_zscore.(currentPeriod);




        % PF correlation and spatial information across quartiles

        for ishuffle = 1:2
            for iseq = 1:2

                currShuffleMethod = shuffles{ishuffle};
                currSeqMetric     = seqMetrics{iseq};

                % PF correlation
                asTuningPFcorr_sub_pooled.(currentPeriod).data.(currShuffleMethod).(currSeqMetric){iSession}   = asTuningPFcorr_sub.(currentPeriod).data.(currShuffleMethod).(currSeqMetric); 
                asTuningPFcorr_sub_pooled.(currentPeriod).ui.(currShuffleMethod).(currSeqMetric){iSession}     = asTuningPFcorr_sub.(currentPeriod).ui.(currShuffleMethod).(currSeqMetric);

                asTuningPFcorr_sub_zscore_pooled.(currentPeriod).(currShuffleMethod).(currSeqMetric){iSession} = asTuningPFcorr_sub_zscore.(currentPeriod).(currShuffleMethod).(currSeqMetric);




                % Spatial information
                asTuningSpatialInfo_sub_pooled.(currentPeriod).data.(currShuffleMethod).(currSeqMetric){iSession} = asTuningSpatialInfo_sub.(currentPeriod).data.(currShuffleMethod).(currSeqMetric);
                asTuningSpatialInfo_sub_pooled.(currentPeriod).ui.(currShuffleMethod).(currSeqMetric){iSession}   = asTuningSpatialInfo_sub.(currentPeriod).ui.(currShuffleMethod).(currSeqMetric);              

                asTuningSpatialInfo_sub_zscore_pooled.(currentPeriod).(currShuffleMethod).(currSeqMetric){iSession} = asTuningSpatialInfo_sub_zscore.(currentPeriod).(currShuffleMethod).(currSeqMetric);

            end
        end
    
    end
    
end




%%

figure;

allCorrValues = [cell2mat(assemblyTuningPFcorr_pooled.PRE.data); cell2mat(assemblyTuningPFcorr_pooled.RUN.data); cell2mat(assemblyTuningPFcorr_pooled.POST.data)];
bins = linspace(min(allCorrValues), max(allCorrValues)+0.01, 15);

colors = {'b'; 'none';'r'};
edgeColors = {'none'; 'k'; 'none'};
for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    h = histc(cell2mat(assemblyTuningPFcorr_pooled.(currPeriod).data(5)), bins); h(end) = [];
    a(iperiod ) = bar(bins(1:end-1)+diff(bins(1:2))/2, h/sum(h), 'Edgecolor', edgeColors{iperiod}, 'facecolor', colors{iperiod}, 'facealpha', 0.3);
    hold on
    
    rr = ylim;
    medianCorr = nanmedian(cell2mat(assemblyTuningPFcorr_pooled.(currPeriod).data(5)));
    
    
    shuffle = cell2mat(assemblyTuningPFcorr_pooled.(currPeriod).ui(5));
    medianCorr_shuffle = nanmedian(shuffle(:));
    
    cl = colors{iperiod}; 
    if iperiod == 2; cl = 'k'; end
 
    line([medianCorr_shuffle medianCorr_shuffle], [0 rr(2)], 'color', cl, 'linewidth', 2, 'linestyle', ':')
    line([medianCorr medianCorr], [0 rr(2)], 'color', cl, 'linewidth', 2)
    
end

xlim([-1 1])
% set(gca,'YColor','none')
% xticklabels([])

rr = legend(a, 'PRE', 'RUN', 'POST');
rr.FontSize = 14;

set(gca, 'box', 'off', 'fontsize', 12)

xlabel({'Assembly tuning - PF correlation'}, 'fontsize', 14)
ylabel('Ratio of units', 'fontsize', 14)
set(gca, 'box', 'off')




%%   Assembly tunings versus BD replay scores     
        
% plot figures


% correlation of assembly tunings with place fields


figure;
set(gcf, 'position', [3000 50 410 900])
colors = {'b'; 'k'; 'r'};
edgeColors = {'k', 'g', 'k'};
panelNumbers = [3 4 5; 7 8 9; 11 12 13];

periodNames = {'PRE'; 'RUN'; 'POST'};


currShuffleMethod = 'unitIDshuffle';
currReplayMethod  = 'replayScore';


for iperiod = 1:3

    currPeriod = periodNames{iperiod};

    subplot(14,1, panelNumbers(iperiod, :))
    set(gca, 'fontsize', 12)
    
    h = violinplot(cell2mat(asTuningPFcorr_sub_pooled.(currPeriod).data.(currShuffleMethod).(currReplayMethod)));

    for ii = 1:4

        h(ii).ViolinColor = colors{iperiod};
        h(ii).EdgeColor   = edgeColors{iperiod};
        h(ii).ShowData    = 0;
    end

    xticklabels({'0-25%'; '25-50%'; '50-75%'; '75-100%'})
    
    if iperiod == 2
        ylabel('assembly tuning - PF correlation', 'fontsize', 18)
    end
    if iperiod == 3; xlabel('Bayesian replay scores', 'fontsize', 14); end
    xtickangle(45)


    ylim([-1 1])
    
    grid on


    hold on
    
    data = cell2mat(asTuningPFcorr_sub_pooled.(currPeriod).data.(currShuffleMethod).(currReplayMethod));
    shuffle = cell2mat(asTuningPFcorr_sub_pooled.(currPeriod).ui.(currShuffleMethod).(currReplayMethod));
    for ii = 1:4
        
        
        temp = shuffle(:, ii, :);
        corrWithPFs_ui = temp(:);
%         validUnits = ~isnan(asTuningPFcorr_sub.(currPeriod)(:, ii));
        dataMedian    = nanmedian(data(:, ii));
        shuffleMedian = nanmedian(corrWithPFs_ui);
        line([ii-0.3 ii+0.3], [shuffleMedian shuffleMedian], 'color', colors{iperiod}, 'linestyle', ':', 'linewidth', 2)

        pval = ranksum(data(:, ii), corrWithPFs_ui, 'tail', 'right');
        
        if pval < 0.001
            text(ii-0.5, -0.5, sprintf('p<0.001'), 'color', [0.1 0.1 0.1], 'fontsize', 8)
        else
            text(ii-0.5, -0.5, sprintf('p = %.3f', pval), 'color', [0.1 0.1 0.1], 'fontsize', 8)
        end
%         text(ii-0.5, dataMedian+0.4, sprintf('n = %d', nPBEs_subset.(currPeriod)(ii)), 'color', 'b', 'fontsize', 8)
    end

    hold off

    tt = ylim;
    text(0.5, tt(2), currPeriod, 'fontsize', 12, 'fontweight', 'bold', 'horizontalAlignment', 'center')

    dim = [0.1 0.8 0.5 .1];
    annotation('textbox',dim,'String', {sessionName; sprintf('%s - %s', currShuffleMethod, currReplayMethod)}, ...
          'FitBoxToText','on', 'Interpreter', 'none');

end

%%

figure; 

periodNames = {'PRE'; 'RUN'; 'POST'};

for ii = 1:3
    
    currPeriod = periodNames{ii};
    ax(ii) = subplot(1,3,ii);

    currAT = assemblyTunings.(currPeriod).data;
    currAT = currAT ./ repmat(max(currAT, [], 2), [1 size(currAT, 2)]);
    
    [nUnits, nPosBins] = size(currAT);
    
    imagesc(2*(1:nPosBins), 1:nUnits, currAT)
    colormap(ax(ii), 'hot')
    xlabel('track position', 'fontsize', 12); 
    if ii == 1 
        ylabel('unit', 'fontsize', 12)
    end
    title([currPeriod sprintf('(n=%d)', nPBEs.(currPeriod))]);
    caxis([0 1])
end


figure; 

spatialTunings = spatialTunings(sortedActiveUnitIndices, :);
spatialTunings = spatialTunings ./ repmat(max(spatialTunings, [], 2), [1 size(currAT, 2)]);

[nUnits, nPosBins] = size(spatialTunings);
    
imagesc(2*(1:nPosBins), 1:nUnits, spatialTunings)
colormap('jet')
xlabel('track position', 'fontsize', 12); 
if ii == 1 
    ylabel('unit', 'fontsize', 12)
end
title([currPeriod sprintf('(n=%d)', nPBEs.(currPeriod))]);
caxis([0 1])



