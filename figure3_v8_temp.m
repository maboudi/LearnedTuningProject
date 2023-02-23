
clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 12 15:17 6];
nSessions    = numel(currSessions);

nShuffles = 1000;

startT = cell(nSessions, 1);
endT   = cell(nSessions, 1);


activeUnits = cell(nSessions, 1);
epochNames = {'pre'; 'run'; 'post'};
dataTypes  = {'data', 'ui'};


% initialize the variables
spikes_pooled         = cell(nSessions, 1);
binCenters            = cell(nSessions, 1);
binCenters_nonOverlap = cell(nSessions, 1);

assemblyTunings_time            = cell(nSessions, 1);
assemblyTunings_time_nonOverlap = cell(nSessions, 1);



% Pearson correlation; between learned tunnings of different units and between learned tunings
% and the correponding place fields of individual units

assemblyTuningPFcorr        = cell(nSessions, 1);
% assemblyTuningPFcorr_pscore = cell(nSessions, 1);

epochsXCorr         = cell(nSessions, 1);
epochsXCorr_zscore  = cell(nSessions, 1);


% KL  diveregenc (the same manner as the Pearson correlation)

assemblyTuningPFKLdiv        = cell(nSessions, 1);
assemblyTuningPFKLdiv_pscore = cell(nSessions, 1);

epochsKLdiv                  = cell(nSessions, 1);
epochsKLdiv_zscore           = cell(nSessions, 1);

load(fullfile(parentDir, 'assemblyTunings_plus_PFunitIDshuffles.mat'), 'learnedTuningPFcorr')


%% calculate matching scores for all sessions

for iSess = 1:nSessions

    
sessionNumber = currSessions(iSess);
sessionName   = rr(sessionNumber+2).name
basePath      = fullfile(parentDir, sessionName);


% % % spikes and behavior data
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'fileInfo', 'spikes_pyr')

behavior = fileInfo.behavior;
behavior.time = behavior.time/3600; % to hours

startT{iSess}.pre  = behavior.time(1,1); endT{iSess}.pre  = behavior.time(2,1); 
startT{iSess}.run  = behavior.time(2,1); endT{iSess}.run  = behavior.time(2,2); 
startT{iSess}.post = behavior.time(2,2); endT{iSess}.post = behavior.time(3,1)+4; % in long recordings only the first 4 hours of post sleep is used in this figures 


spikes = spikes_pyr;

nUnits   = numel(spikes);
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);



% bidirectional spatial tuning

spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;    
end

spikes_pooled{iSess} = spikes;



% % % SWR LEARNED TUNINGS


% % learned tunings calculated based on the entire epochs
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), ...
    'activeUnits', ...
    'assemblyTuningPFcorr');%, ...
%     'assemblyTuningPFKLdiv');

activeUnits{iSess} = s.activeUnits;

okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


% assemblyTuningPFcorr{iSess}  = s.assemblyTuningPFcorr;


for iEpoch = 1:3
    for idt = 1:2
    
        assemblyTuningPFcorr{iSess}.(epochNames{iEpoch}).(dataTypes{idt}) = learnedTuningPFcorr.PBEs.(epochNames{iEpoch}).(dataTypes{idt}){iSess};

    end
end





% assemblyTuningPFKLdiv{iSess} = s.assemblyTuningPFKLdiv;

% for iepoch = 1:3
%     currEpoch = epochNames{iepoch};
%     
%     allPFcorr_prctiles = calPrctile(s.assemblyTuningPFcorr.(currEpoch));
%     assemblyTuningPFcorr_pscore{iSess}.(currEpoch) = allPFcorr_prctiles;
%     
% %     allKLdiv_prctiles = calPrctile(s.assemblyTuningPFKLdiv.(currEpoch));
% %     assemblyTuningPFKLdiv_pscore{iSess}.(currEpoch) = allKLdiv_prctiles;
% 
% end




% % learned tunings in 15-minutes time windows

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time_Jan2022.mat']), ...
        'binCenters', 'assemblyTunings_time');
    
binCenters{iSess}           = s.binCenters;
assemblyTunings_time{iSess} = s.assemblyTunings_time;

for iEpoch = 1:numel(epochNames)
    binCenters{iSess}.(epochNames{iEpoch}) = binCenters{iSess}.(epochNames{iEpoch})/3600;
end



% truncate the long recording sessions to the first 4 hours of post sleep

idx = find(binCenters{iSess}.post < endT{iSess}.post);

binCenters{iSess}.post = binCenters{iSess}.post(idx);
assemblyTunings_time{iSess}.post.data = assemblyTunings_time{iSess}.post.data(:, :, idx);




% remove learned tunnigs in overlapping time windows
% we used 15minutes-long sliding windows with 5minutes step size, so we
% need to keep every third time windows
    
for iEpoch = 1:numel(epochNames)
    epochName = epochNames{iEpoch};
    
    binCenters_nonOverlap{iSess}.(epochName) = binCenters{iSess}.(epochName)(1:3:end);
    assemblyTunings_time_nonOverlap{iSess}.(epochName).data = assemblyTunings_time{iSess}.(epochName).data(:, :, 1:3:end);
end




%% calculate a cross-correlation matrix for the 15minues time window learned tunings
% removed the part related to KL divergence

epochsXCorr{iSess}.data   = nan(3, 3, nUnits);
epochsXCorr{iSess}.ui     = nan(3, 3, nShuffles);
epochsXCorr_zscore{iSess} = nan(3, 3, nUnits); 
   

for iEpoch = 1:3
    
    currData1 = assemblyTunings_time_nonOverlap{iSess}.(epochNames{iEpoch}).data;

    for jEpoch = iEpoch:3
        currData2 = assemblyTunings_time_nonOverlap{iSess}.(epochNames{jEpoch}).data;
        
        %%% temporary code to check if the pre LT stability biases its
        %%% similarity with post LT
        
        
%         if jEpoch ~= iEpoch
%             currData2_2 = nan(size(currData2));
% 
% 
%             randOkUnits             = okUnits(randperm(numel(okUnits)));
%             currData2_2(okUnits, :) = currData2(randOkUnits, :);
% 
%             currData2 = currData2_2;
%         end
        
        %%%

        for iUnit = 1:nUnits
            
%             otherUnits = setdiff(1:nUnits, iUnit);
                        
            % % data
            
            % the real tunings
            tuning1 = squeeze(currData1(iUnit, :, :));
            tuning2 = squeeze(currData2(iUnit, :, :));

            
            % calculate the matching scores: pearson correlation and KL
            % divergence
            corrMat = corr(tuning1, tuning2, 'type', 'pearson');
            
            if jEpoch == iEpoch
                
                offDiagonalIdx = true(size(corrMat));
                for ii = 1:size(corrMat, 1)
                    offDiagonalIdx(ii, ii) = false;
                end
                corrMat = corrMat(offDiagonalIdx);
            end
            
            epochsXCorr{iSess}.data(iEpoch, jEpoch, iUnit) = nanmedian(corrMat(:));

        end
        
        
        % % shuffle 
        for inst = 1:nShuffles

            % generate shuffle tunings
            nTimeBins1 = size(currData1, 3);
            tuning1    = zeros(nPosBins, nTimeBins1);
            for iTimeBin = 1:nTimeBins1
                tuning1(:, iTimeBin) = squeeze(currData1(okUnits(randperm(numel(okUnits), 1)), :, iTimeBin));
            end
            

            nTimeBins2 = size(currData2, 3);
            tuning2 = zeros(nPosBins, nTimeBins2);
            for iTimeBin = 1:nTimeBins2 
                tuning2(:, iTimeBin) = squeeze(currData2(okUnits(randperm(numel(okUnits), 1)), :, iTimeBin));
            end


            % calculate the matching scores
            corrMat = corr(tuning1, tuning2, 'type', 'pearson');

            if jEpoch == iEpoch
                offDiagonalIdx = true(size(corrMat));
                for ii = 1:size(corrMat, 1)
                    offDiagonalIdx(ii, ii) = false;
                end
                corrMat = corrMat(offDiagonalIdx);
            end

            epochsXCorr{iSess}.ui(iEpoch, jEpoch, inst) = nanmedian(corrMat(:));

        end
    end    
end

epochsXCorr_zscore{iSess} = ((epochsXCorr{iSess}.data) - nanmean(epochsXCorr{iSess}.ui, 3)) ./ nanstd(epochsXCorr{iSess}.ui, [], 3);

end


% %% plot the 15-minutes time windows learned tunings and matching scores for all units in a session
% 
% 
% plotLearnedTunings_allUnits




%% plot the learning tunings across time for example units in an example session

iSess = 2;
% selectedUnits = [72 117 79 116 63 90]; % AchillesMaze1
selectedUnits = [92 40 65 74 96 61]; % AchillesMaze2
% selectedUnits = [33 46 85 90 123 150]; % rat U


% selectedUnits = [111 118 119 130 135]; % Achilles_10252013
% selectedUnits = [15 60]; %Achilles_112520113
% selectedUnits = [9 34]; %Cicero
% selectedUnits = [7 21 26 29 45];% Gatsby
% selectedUnits = [8];


spikes = spikes_pooled{iSess};



plotwidth  = 360;
plotheight = 235;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

nor = 3;
noc = 2;

leftmargin = 30;  rightmargin = 30;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 8;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);

nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

for iUnit = 1:numel(selectedUnits)
    currUnit = selectedUnits(iUnit); 
    
    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on

    currSpatialTuning = spikes(currUnit).spatialTuning_smoothed.uni;

    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time{iSess}.(currEpoch).data(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch}./repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

        imagesc(binCenters{iSess}.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); colormap('jet')

    end
    
    plot(0.1*currSpatialTuning + startT{iSess}.run , 1:nPosBins, 'linewidth', 1.5, 'color', [0 0 0 0.6], 'DisplayName', 'track spatial tuning')


    yl = ylim(ax(iUnit));

    
    if iUnit == 4
        line(0.2*[0 2] + startT{iSess}.run, [yl(2)-10 yl(2)-10], 'linewidth', 1.5, 'color', [0 0 0 0.6])
        text(startT{iSess}.run, yl(2)+20, '2 Hz', 'fontsize', fontsize, 'color', 'k')
    end
    
    ylim([0 nPosBins]);
    yl = ylim;
    text(0, yl(2)+12, ['unit ' num2str(iUnit)], 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'left')
    
    xlim([0 binCenters{iSess}.post(end)])
    
    set(ax(iUnit), 'YTick', [1 nPosBins], 'YTickLabel', [0 1], 'TickDir', 'out','TickLength',[0.01, 0.01])
    
    if ismember(iUnit, [1 4])
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           currBinCenters = binCenters{iSess}.(currEpoch);
           text(median(currBinCenters), yl(2)+20, currEpoch, 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'center')
       end
    end
        
    
    if iUnit == 2
        ylabel('position(normalized)', 'fontsize', fontsize)
    end
    
    if ismember(iUnit, [3 6])
        xlabel('time (hour)')
    end
    
    if iUnit >= 4
        set(ax(iUnit), 'YTickLabel', {})
    end
    
    if ~ismember(iUnit, [3 6])
       set(ax(iUnit), 'XTickLabel', {}) 
    end
    
    hold off
    
end



%% plot cross-correlation matrix of learned tunings for the entire session


plotwidth  = 375;
plotheight = 175;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nor = 1;
noc = 6;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


cnctBinCenters = [binCenters{iSess}.pre binCenters{iSess}.run binCenters{iSess}.post];
tickLabels = 0:2:cnctBinCenters(end);

for iUnit = 1:numel(selectedUnits)
    
    
    currUnit = selectedUnits(iUnit);
    
    
    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time{iSess}.(currEpoch).data(currUnit, :, :, 1));

    end
    
    tmpConcat = cell2mat(learnedTuning');

    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    
    
    for iEpoch = 1:3
        for jEpoch = 1:3
            
            corrMat = corr(learnedTuning{jEpoch}, learnedTuning{iEpoch});
            corrMat(isnan(corrMat)) = 0;
            imagesc(binCenters{iSess}.(epochNames{iEpoch}), binCenters{iSess}.(epochNames{jEpoch}), corrMat)
 
        end
    end
    caxis([-1 1])
    colormap('jet')
    

    ylim([0 cnctBinCenters(end)])
    xlim([0 cnctBinCenters(end)])
    
    yl = ylim;
    text(0, yl(2)+0.5, sprintf('unit %d', iUnit), 'fontsize', fontsize, 'HorizontalAlignment', 'left')
    
    if iUnit >= 2
        set(gca, 'YTickLabel', {})
    end
    
    if iUnit == 1
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           h = text(-4, median(binCenters{iSess}.(currEpoch)), currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           h = text(-2, median(cnctBinCenters), 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           
           text(median(binCenters{iSess}.(currEpoch)), -4 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           
          text(median(cnctBinCenters), -2 , 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 

       end
    end
    
    
    xticks(tickLabels)
    yticks(tickLabels)
    axis square
    
    if iUnit == 6
        h = colorbar;
        h.Label.String = 'correlation';
    end
    set(gca, 'box', 'off', 'fontsize', fontsize, 'position', sub_pos{iUnit}, 'TickDir', 'out','TickLength',[0.01, 0.01])
   
    
end




%% plot the distribution of mean cross-correlations for the example units


% for the example units (distribution of correlation for real data and ui)

plotwidth  = 330;
plotheight = 150;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nor = 1;
noc = 6;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


for iUnit = 1:numel(selectedUnits)
    
    currUnit = selectedUnits(iUnit);
    
    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on
    
    
    tt = 0;
    for iEpoch = [1 3] % pre and post within-epoch consistency of learned tunings over time
        tt = tt + 1; 
        % data
        data = epochsXCorr{iSess}.data(iEpoch, iEpoch, currUnit);
        
%         [xData, yData] = myViolin(data);
% 
%         patch(tt-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        scatter(tt-0.15, data, 5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')

        % ui
        ui = squeeze(epochsXCorr{iSess}.ui(iEpoch, iEpoch, :));
        
        
        try 
        [xData, yData] = myViolin(ui);
        patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
        scatter(tt+0.15, nanmedian(ui), 2, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')

        catch
            continue
        end
        
    end
    
    
    tt = tt+1;
    
    % similarity of learned tunings between pre and post
    
    % data
    data = epochsXCorr{iSess}.data(1, 3, currUnit, 1);
        
%     [xData, yData] = myViolin(data);
% 
%     h1 = patch(tt-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    h1 = scatter(tt-0.15, data, 5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none'); 
    
    
    % ui
    ui = squeeze(epochsXCorr{iSess}.ui(1, 3, :));
    
    try
        [xData, yData] = myViolin(ui);    
        h2 = patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        scatter(tt+0.15, nanmedian(ui), 2, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none');
    catch
        continue
    end
    
    
    
    if iUnit == 1
        ylabel({'learned tuning'; 'stability index'}, 'fontsize', fontsize)
    end
    
    if iUnit == 6
       legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')
    end
    
    xlim([0.5 3.5])
%     ylim([-1 1])
%     grid on
    
    if iUnit >= 2
        set(gca, 'YTickLabel', {})
    end
    
    set(ax(iUnit), 'fontsize', fontsize, 'position', sub_pos{iUnit}, 'linewidth', 1, 'box', 'off', 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
    xtickangle(ax(iUnit), 45)
    set(gca, 'ytick', -1:0.5:1)
    
    aa = gca;
    aa.YGrid = 'on';
    
end

linkaxes(ax, 'y')
yl = ylim;

ylim([yl(1)-0.1 yl(2)+0.1])




%% plot the distribution of cross-correlation for an example session (panel D)

iSess = 2;

plotwidth  = 120;
plotheight = 120;
fontsize   = 6;

okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

med     = nan(3,1);
p_value = nan(3, 1);
z_value = nan(3, 1);

tt = 1;
for iEpoch = [1 3]
    
   
    currMeds = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits, 1));
   
    scatter(tt+0.1*randn(numel(currMeds), 1), currMeds, 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    
    [xData, yData] = myViolin(currMeds);
    patch(tt+5*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    
 
    med(tt) = nanmedian(currMeds);


    [p_value(tt), ~, stats] = signrank(currMeds, 0, 'tail', 'right');
    z_value(tt) = stats.zval;

    line([tt-0.2 tt+0.2], [med(tt) med(tt)], 'linewidth', 1, 'color', 'k')
    
    
    tt = tt + 1;
end

currMeds = squeeze(epochsXCorr_zscore{iSess}(1, 3, okUnits, 1));

scatter(tt+0.1*randn(numel(currMeds), 1), currMeds, 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    
[xData, yData] = myViolin(currMeds);
patch(tt+5*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)

med(tt) = nanmedian(currMeds);

[p_value(tt), ~, stats] = signrank(currMeds, 0, 'tail', 'right');
z_value(tt) = stats.zval;


line([tt-0.2 tt+0.2], [med(tt) med(tt)], 'linewidth', 1, 'color', 'k')   


xlim([0.5 3.5])

% set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:2, 'xticklabel', {'pre'; 'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])

xtickangle(gca, 45)
ylabel('LT stability/similarity (z)', 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';




%% plot the distribution of cross-corrrelation for all sessions

plotwidth  = 170;
plotheight = 120;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


hold on

nSessions = size(epochsXCorr_zscore, 1);
med_pooled = nan(3,1);
p_value    = nan(3,1);
z_value    = nan(3,1);

med= [];

tt = 1;
for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};

    currData.(currEpoch) = cell(nSessions, 1);
    for iSess = 1:nSessions
        curr_x = tt+(iSess-ceil(nSessions/2))*0.85/nSessions;
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        currData.(currEpoch){iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits)); % for plotting the stability
        
%         if iEpoch == 1 % for plotting the matching between different LTs in different epochs
%             currData{iSess} = squeeze(epochsXCorr_zscore{iSess}(1, 3, okUnits)); 
%         elseif iEpoch == 3
%             currData{iSess} = squeeze(epochsXCorr_zscore{iSess}(2, 3, okUnits));
%         end

        med.(currEpoch)(iSess) = nanmedian(currData.(currEpoch){iSess});
        lq  = nanprctile(currData.(currEpoch){iSess}, 25);
        hq  = nanprctile(currData.(currEpoch){iSess}, 75);
        data_min = min(currData.(currEpoch){iSess});
        data_max = max(currData.(currEpoch){iSess});


        scatter(curr_x, med.(currEpoch)(iSess), 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

        % whiskers
        line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 0 0 0.8])
        line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 0 0 0.8])

        % box
        patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    end
    
    pooledData.(currEpoch) = cell2mat(currData.(currEpoch));
    
    [p_value(tt), ~, stats] = signrank(pooledData.(currEpoch), 0, 'tail', 'right');
    z_value(tt) = stats.zval;


    med_pooled(tt) = nanmedian(pooledData.(currEpoch));
    lq  = nanprctile(pooledData.(currEpoch), 25);
    hq  = nanprctile(pooledData.(currEpoch), 75);
    data_min = min(pooledData.(currEpoch));
    data_max = max(pooledData.(currEpoch));

    curr_x = tt + (iSess+1-ceil(nSessions/2))*0.85/nSessions;

    scatter(curr_x, med_pooled(tt), 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

    % whiskers
    line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 1 0 0.8])
    line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 1 0 0.8])

    % box
    patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    
    tt = tt + 1;
end


currData.pre_w_post = cell(nSessions, 1);
for iSess = 1:nSessions
    curr_x = tt+(iSess-ceil(nSessions/2))*0.85/nSessions;
    okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);

    currData.pre_w_post{iSess} = squeeze(epochsXCorr_zscore{iSess}(1, 3, okUnits)); % for plotting the stability

    med.pre_w_post(iSess) = nanmedian(currData.pre_w_post{iSess});
    lq  = nanprctile(currData.pre_w_post{iSess}, 25);
    hq  = nanprctile(currData.pre_w_post{iSess}, 75);
    data_min = min(currData.pre_w_post{iSess});
    data_max = max(currData.pre_w_post{iSess});


    scatter(curr_x, med.pre_w_post(iSess), 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

    % whiskers
    line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 0 0 0.8])
    line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 0 0 0.8])

    % box
    patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
end

pooledData.pre_w_post = cell2mat(currData.pre_w_post);

[p_value(tt), ~, stats] = signrank(pooledData.pre_w_post, 0, 'tail', 'right');
z_value(tt) = stats.zval;


med_pooled(tt) = nanmedian(pooledData.pre_w_post);
lq  = nanprctile(pooledData.pre_w_post, 25);
hq  = nanprctile(pooledData.pre_w_post, 75);
data_min = min(pooledData.pre_w_post);
data_max = max(pooledData.pre_w_post);

curr_x = tt+(iSess+1-ceil(nSessions/2))*0.85/nSessions;

scatter(curr_x, med_pooled(tt), 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

% whiskers
line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 1 0 0.8])
line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 1 0 0.8])

% box
patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8)





xlim([0.5 3.5])

% set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:2, 'xticklabel', {'pre';'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])


xtickangle(gca, 45)
% ylabel({'learned tuning'; 'stability index(z)'}, 'fontsize', fontsize)
ylabel('LT stability/similarity (z)', 'fontsize', fontsize)

ax = gca;
ax.YGrid = 'on';



%% (added later) violin plot of pooled data

iSess = 2;

plotwidth  = 120;
plotheight = 120;
fontsize   = 6;



f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

med_pooled     = nan(3,1);
p_value = nan(3, 1);
z_value = nan(3, 1);

tt = 1;
for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
   
    currMeds = pooledData.(currEpoch);
   
%     scatter(tt+0.1*randn(numel(currMeds), 1), currMeds, 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    indivSessionMedian = med.(currEpoch);
    
    dots = indivSessionMedian(setdiff(1:nSessions, 2));
    scatter(tt+0.1*randn(numel(dots), 1), dots, 2, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    scatter(tt+0.1*randn(1), indivSessionMedian(2), 2, 'MarkerFaceColor', 'g', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);


    [xData, yData] = myViolin(currMeds);
    patch(tt+4*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    
 
    med_pooled(tt) = nanmedian(currMeds);

% 
%     [p_value(tt), ~, stats] = signrank(currMeds, 0, 'tail', 'right');
%     z_value(tt) = stats.zval;
%     
    [p_value(tt), ~, stats] = signrank(indivSessionMedian, 0, 'tail', 'right');
%     z_value(tt) = stats.zval;


    line([tt-0.5 tt+0.5], [med_pooled(tt) med_pooled(tt)], 'linewidth', 1, 'color', 'k')
    
    
    tt = tt + 1;
end

currMeds = pooledData.pre_w_post;

% scatter(tt+0.1*randn(numel(currMeds), 1), currMeds, 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);

indivSessionMedian = med.pre_w_post;

dots = indivSessionMedian(setdiff(1:nSessions, 2));
scatter(tt+0.1*randn(numel(dots), 1), dots, 2, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
scatter(tt+0.1*randn(1), indivSessionMedian(2), 2, 'MarkerFaceColor', 'g', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
   
    
[xData, yData] = myViolin(currMeds);
patch(tt+4*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)

med_pooled(tt) = nanmedian(currMeds);

[p_value(tt), ~, stats] = signrank(indivSessionMedian, 0, 'tail', 'right');% currMeds
% z_value(tt) = stats.zval;


line([tt-0.5 tt+0.5], [med_pooled(tt) med_pooled(tt)], 'linewidth', 1, 'color', 'k')   


xlim([0.5 3.5])

% set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:2, 'xticklabel', {'pre'; 'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])

xtickangle(gca, 45)
ylabel('LT stability/similarity (z)', 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';



% comparing between the different epochs/conditions



[pp, ~, stats] = signrank(pooledData.pre, 0, 'tail', 'right');
zz = stats.zval;





%% correlation between LT consistency during an epoch and its correlation with PFs


plotwidth  = 90;
plotheight = 130;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nSessions = size(epochsXCorr_zscore, 1);
thresh_z = 2;

clear ax med LTstability_pooled PFcorrs_pooled


for iEpoch = 1:3
    
    hold on  
    
    LTstability  = cell(nSessions, 1); % LT consistency
    
%     PFcorrs_zscore = cell(nSessions, 1);
%     KLdivs_zscore  = cell(nSessions, 1);

    
    PFcorrs.data = cell(nSessions, 1);
    PFcorrs.ui   = cell(nSessions, 1);

    KLdivs  = cell(nSessions, 1);
        
    pre_run_matchingIdx  = cell(nSessions, 1);
    pre_post_matchingIdx = cell(nSessions, 1);
    run_post_matchingIdx = cell(nSessions, 1);
    unitIDs              = cell(nSessions, 1);
    
    for iSess = 1:nSessions 
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        unitIDs{iSess} = [iSess*ones(numel(okUnits), 1) okUnits];

        
        LTstability{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));

        PFcorrs.data{iSess} = assemblyTuningPFcorr{iSess}.(epochNames{iEpoch}).data(okUnits);
        PFcorrs.ui{iSess}   = assemblyTuningPFcorr{iSess}.(epochNames{iEpoch}).ui(okUnits, :);

%         KLdivs{iSess}  = assemblyTuningPFKLdiv{iSess}.(epochNames{iEpoch}).data(okUnits);
        
        
%         PFcorrs_zscore{iSess} = assemblyTuningPFcorr_zscore{iSess}.(epochNames{iEpoch})(okUnits);
%         KLdivs_zscore{iSess}  = assemblyTuningPFKLdiv_zscore{iSess}.(epochNames{iEpoch})(okUnits);

        pre_run_matchingIdx{iSess}  = squeeze(epochsXCorr_zscore{iSess}(1, 2, okUnits));
        pre_post_matchingIdx{iSess} = squeeze(epochsXCorr_zscore{iSess}(1, 3, okUnits));
        run_post_matchingIdx{iSess} = squeeze(epochsXCorr_zscore{iSess}(2, 3, okUnits));
       
        
    end

    if iEpoch == 1
        unitIDs_pooled = cell2mat(unitIDs);
    end
    
    LTstability_pooled.(epochNames{iEpoch}) = cell2mat(LTstability);

%     KLdivs_zscore_pooled.(epochNames{iEpoch})  = cell2mat(KLdivs_zscore);
%     PFcorrs_zscore_pooled.(epochNames{iEpoch}) = cell2mat(PFcorrs_zscore);
    
%     KLdivs_pooled.(epochNames{iEpoch})  = cell2mat(KLdivs);
    PFcorrs_pooled.(epochNames{iEpoch}).data = cell2mat(PFcorrs.data);
    PFcorrs_pooled.(epochNames{iEpoch}).ui   = cell2mat(PFcorrs.ui);
    
end

pre_run_matchingIdx_pooled  = cell2mat(pre_run_matchingIdx);
pre_post_matchingIdx_pooled = cell2mat(pre_post_matchingIdx);
run_post_matchingIdx_pooled = cell2mat(run_post_matchingIdx);




tt = 1;
for iEpoch = [3]
    
    currEpoch = epochNames{iEpoch};
    
    currParam = PFcorrs_pooled.(currEpoch);
%     currParam = pre_post_matchingIdx_pooled;
    
    stableOnes    = LTstability_pooled.(currEpoch) > thresh_z;
    nonStableOnes = LTstability_pooled.(currEpoch) <= thresh_z;
    
    
%     pooledData = currParam; %[currParam.pre; currParam.run; currParam.post];
    
%     bottomLim = prctile(pooledData, 0);
%     topLim    = prctile(pooledData, 100);
    
    
    stableData = currParam.data(stableOnes); %.(epochNames{iEpoch})
    % check if the stable data has a median above corresp. shuffle (added later)
    
    med_data = nanmedian(currParam.data(stableOnes));
    med_ui   = nanmedian(currParam.ui(stableOnes, :));

    pvalue_stable = numel(find(med_ui >= med_data))/size(med_ui, 2);

    %

%     stableData = stableData(stableData >= bottomLim & stableData <= topLim);
    [xData, yData] = myViolin_oneSided(stableData, 0.01);
    h1 = patch(tt+5*xData, yData, [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    lq  = prctile(stableData, 25);
    hq  = prctile(stableData, 75);
    med = median(stableData);
    
    line([tt tt+0.5], [lq lq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
    line([tt tt+0.5], [hq hq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
    line([tt tt+0.5], [med med], 'linestyle', '-', 'color', 'w', 'linewidth', 1)


    nonStableData = currParam.data(nonStableOnes); %.(epochNames{iEpoch})
    % check if the nonStable data had a median above corresp. shuffle (added later)

    med_data = nanmedian(currParam.data(nonStableOnes));
    med_ui   = nanmedian(currParam.ui(nonStableOnes, :));

    pvalue_nonStable = numel(find(med_ui >= med_data))/size(med_ui, 2); 


    %


%         nonStableData = nonStableData(nonStableData >= bottomLim & nonStableData <= topLim);
    [xData, yData] = myViolin_oneSided(nonStableData, 0.01); 
    h2 = patch(tt-5*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    lq  = prctile(nonStableData, 25);
    hq  = prctile(nonStableData, 75);
    med = median(nonStableData);

    line([tt-0.5 tt], [lq lq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
    line([tt-0.5 tt], [hq hq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
    line([tt-0.5 tt], [med med], 'linestyle', '-', 'color', 'w', 'linewidth', 1)

    

    [pval, ~,stats] = ranksum(stableData, nonStableData, 'tail', 'right');
    z_value = stats.zval;
%     med.(epochNames{iEpoch}) = [nanmedian(currParam.(epochNames{iEpoch})(nonStableOnes)) nanmedian(currParam.(epochNames{iEpoch})(stableOnes))];
    
    tt = tt+1;
end

ax = gca;
ax.YGrid = 'on';

xlim([0.5 2.5])


% legend([h1, h2], 'stable LTs', 'nonstable LTs', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')

% xlabel('epoch', 'fontsize', fontsize)
ylabel('PF-LT Pearson correlation coeff. (z)', 'fontsize', fontsize)
% ylabel('PF-LT KL divergence (z)', 'fontsize', fontsize)


yl = ylim;

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1, 'xticklabel', {'pre'}, 'ytick', yl(1):1:yl(2), 'TickDir', 'out','TickLength',[0.01, 0.01] ) 


    

%% plot pre_post_matchingIdx versus post LT-PF correlation for units with
% stable LTs during the post


plotwidth  = 120;
plotheight = 130;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

idx = LTstability_pooled.post > 8;
aa = pre_post_matchingIdx_pooled(idx);
bb = PFcorrs_pooled.post(idx);
scatter(bb, aa, 2, 'filled', 'markerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)


X = [ones(size(aa)) bb];
[b, ~,~,~, stats] = regress(aa, X);
pval = stats(3);
r2   = stats(1);

aaCalc = b(1)+b(2)*bb;
plot(bb, aaCalc, 'color', [0 0 0 0.5], 'linewidth', 1)
text(bb(end), aaCalc(end), {sprintf('R2=%.2f', r2); ''; sprintf('p=%.1e', pval)}, 'fontsize', fontsize)

ylabel('pre-post LT similarity (z)', 'fontsize', fontsize)
xlabel('post LT-PF similarity', 'fontsize', fontsize)

set(gca, 'box', 'off', 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01], 'linewidth', 1)


%% how the r2 value calculated above changes as a function of Post stability

stabilityThreshs = [-inf 2 4 6 8 10 12];

n = numel(stabilityThreshs);
pval = nan(n, 1);
r2   = nan(n, 1);

for is = 1:n
    
    idx = LTstability_pooled.post > stabilityThreshs(is);
    aa = pre_post_matchingIdx_pooled(idx);
    bb = PFcorrs_pooled.post(idx);
    
    X = [ones(size(bb)) aa];
    [b, ~,~,~, stats] = regress(bb, X);
    pval(is) = stats(3);
    r2(is)   = stats(1);
    
end

plotwidth  = 80;
plotheight = 70;
fontsize   = 4;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on
plot([0 stabilityThreshs(2:end)], r2, 'color', [0 0 0 0.5], 'marker', 'square', 'markersize', 2, 'linewidth', 0.5)

x = [0 stabilityThreshs(2:end)];
for is = 1:n
    text(x(is), r2(is), sprintf('p=%.1e', pval(is)), 'fontsize', 4)
end

xlabel('stability threshold', 'fontsize', fontsize)
ylabel('R2', 'fontsize', fontsize)

set(gca, 'xtick', [0 stabilityThreshs(2:end)], 'TickDir', 'out','TickLength',[0.01, 0.01], 'linewidth', 0.5, 'box', 'off', 'fontsize', fontsize)

%%


plotwidth  = 100;
plotheight = 110;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on


cc = LTstability_pooled.post;


aa = LTstability_pooled.post;
bb = LTstability_pooled.pre;



% scatter(bb, aa, 2, 'filled', 'markerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
scatter(bb, aa, 1, 'k', 'filled', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
% colormap('jet')
% h = colorbar('orientation', 'Horizontal', 'location', 'north');
% h.Label.String = 'post LT stability';
% caxis(old_caxis)


X = [ones(size(aa)) bb];
[b, ~,~,~, stats] = regress(aa, X);
pval = stats(3);
r2   = stats(1);

aaCalc = b(1)+b(2)*bb;
h1 = plot(bb, aaCalc, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(bb)-5, nanmedian(aaCalc), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)


minn = min([aa; bb]);
maxx = max([aa; bb]);

h2 = line([minn maxx], [minn maxx], 'color', [0 0 0 0.5], 'linewidth', 1, 'linestyle', ':');



xlabel('pre LT stability (z)', 'fontsize', fontsize)
ylabel('post LT stability (z)', 'fontsize', fontsize)


xlim([-12 40])
ylim([-12 40])
set(gca, 'XTick', -10:10:40, 'YTick', -10:10:40, 'box', 'off', 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01], 'linewidth', 1)

legend([h1, h2], {'regression line'; 'identity line'}, 'box', 'off', 'fontsize', 4, 'location', 'northout')

axis square


%% pre_post similarity versus LT stability in pre or post

plotwidth  = 120;
plotheight = 130;
fontsize   = 12;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

variable = LTstability_pooled.pre;
binRanges = floor(min(variable)):7:ceil(max(variable));
% binRanges(end) = ceil(max(variable));

binCenters = binRanges(1:end-1)+1;

[binCounts, idx] = histc(variable, binRanges);
nBins = numel(binRanges)-1;

xticklabels = cell(nBins, 1);
for ii = 1:nBins
    xticklabels{ii} = sprintf('%d-%d', binRanges(ii), binRanges(ii+1));
end



boxplot(LTstability_pooled.post, idx, 'colors', [0.3 0.3 0.3], 'whisker', inf)

set(gca, 'XTickLabel', xticklabels)

xlabel('pre LT stability (z)', 'fontsize', fontsize)
ylabel('post LT stability (z)', 'fontsize', fontsize)

set(gca, 'box', 'off','linewidth', 1, 'fontsize', fontsize, ...
             'TickDir', 'out','TickLength',[0.01, 0.01])

xtickangle(gca, 45)
         
         


% 
% nBins = numel(binRanges)-1;
% 
% for ii = 1:nBins
%     
%      = pre_post_matchingIdx_pooled(idx == ii);








%% a small analysis on correlation between LT-PF corrs and LT-PF KLdivs


KLdivs_all = [KLdivs_zscore_pooled.pre  KLdivs_zscore_pooled.post]; 
PFcorrs_zscore_all = [PFcorrs_zscore_pooled.pre  PFcorrs_zscore_pooled.post];

idx =  PFcorrs_zscore_all > prctile(PFcorrs_zscore_all, 0) & PFcorrs_zscore_all < prctile(PFcorrs_zscore_all, 100); % KLdivs_all > prctile(KLdivs_all, 1) & KLdivs_all < prctile(KLdivs_all, 99) &

KLdivs_all = KLdivs_all(idx);
PFcorrs_zscore_all = PFcorrs_zscore_all(idx);



plotwidth  = 150;
plotheight = 150;
fontsize   = 9;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


scatter(KLdivs_all, PFcorrs_zscore_all, 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.4)

xlabel('LT-PF KL divergence', 'fontsize', fontsize)
ylabel('LT-PF Pearson correlation', 'fontsize', fontsize)

% title('correlation between PF matching scores', 'fontsize', 6, 'fontweight', 'normal')

[coeff, pval] = corr(KLdivs_all, PFcorrs_zscore_all, 'type', 'pearson');

xl = xlim;
yl = ylim;

text(median(xl)-0.2*range(xl), yl(1)+ 0.2*range(yl), {sprintf('Pearson corr coeff. = %.2f', coeff); sprintf('pval = %.2f', pval)}, 'fontsize', 6, 'color', [0.4 0.4 0.4])

set(gca, 'box', 'off', 'linewidth', 1, 'xtick', [xl(1):1:xl(2)], 'ytick', [yl(1):1:yl(2)], 'TickDir', 'out','TickLength',[0.01, 0.01] )

grid on




%% correlation between LT consistency during an epoch and its correlation with PFs

plotwidth  = 200;
plotheight = 250;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



nor = 2;
noc = 1;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 0;
gapr = 5;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


nSessions = size(epochsXCorr_zscore, 1);


clear ax 

tt = 1;
for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
    
    ax(tt) = axes('position',sub_pos{tt},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on
        
    LTstability     = cell(nSessions, 1); % LT consistency
    PFcorrs_zscore         = cell(nSessions, 1);
    PFcorrs_zscored = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTstability{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
        PFcorrs{iSess} = assemblyTuningPFcorr_pooled{iSess}.(currEpoch)(okUnits);
        PFcorrs_zscored{iSess} = assemblyTuningPFcorr_zscore_pooled{iSess}.(currEpoch)(okUnits);
    end
    
    LTstability_pooled.(currEpoch)     = cell2mat(LTstability);
    PFcorrs_pooled.(currEpoch)         = cell2mat(PFcorrs);
    PFcorrs_zscored_pooled.(currEpoch) = cell2mat(PFcorrs_zscored);
    
    
    consisUnits = find(abs(LTstability_pooled.(currEpoch)) > 2);
    
    
%     allCorr = PFcorrs_pooled.(currEpoch)(consisUnits);
    signIdx = abs(PFcorrs_zscored_pooled.(currEpoch)(consisUnits)) > 1 & PFcorrs_zscored_pooled.(currEpoch)(consisUnits).*PFcorrs_pooled.(currEpoch)(consisUnits) > 0; 
    
    
    binEdges   = -1:0.05:1;
    binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;
    
    allCounts    = histc(PFcorrs_pooled.(currEpoch)(consisUnits), binEdges); allCounts(end) = [];
    nonSigCounts = histc(PFcorrs_pooled.(currEpoch)(consisUnits(~signIdx)), binEdges); nonSigCounts(end) = [];
    

    bar(binCenters, allCounts, 'FaceColor', 'k', 'FaceAlpha', 0.8, 'EdgeColor', 'k')
    bar(binCenters, nonSigCounts, 'FaceColor', 'w', 'EdgeColor', 'k')
    
    
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 
    
    tt = tt + 1;

end

% 
% xlim([0.5 2.5])
% 
% 
% xlabel('learned tuning consistency', 'fontsize', fontsize)
% ylabel('PF correlation', 'fontsize', fontsize)


    
%%

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


function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end



function [xData, yData] = myViolin(data)

binEdges = linspace(min(data)-eps, max(data)+eps, 50);

count = histc(data, binEdges); 
count(end) = [];
binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;


count = count/sum(count);
count = count';

gw = gausswindow(2,6);
count = conv(count, gw, 'same');

xData = [count -fliplr(count)];
yData = [binCenters fliplr(binCenters)];

end


function [xData, yData] = myViolin_oneSided(data, binSize)

% data = data(data > prctile(data, 5) & data < prctile(data, 97.5));

binEdges = (min(data)-eps):binSize:(max(data)+eps);
count = histc(data, binEdges); 
count(end) = [];
binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;


count = count/sum(count);
count = count';

gw = gausswindow(6,18);
count = conv(count, gw, 'same');

xData = [count zeros(1, numel(count))];
yData = [binCenters fliplr(binCenters)];

end