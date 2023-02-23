clear; clc; close all

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

sub = dir(parentDir);
nSessions    = numel(sub) - 2; % subtracting two becuase of '.' and '..'
% sessionNames = cell(nSessions, 1);


epochNames = {'pre'; 'run'; 'post'};


replayScoreMethods = {'rt_ui'; 'wc_ui'; 'wc_ts'};


replayScoreMethods_fullName = {'radon Integral - unit ID shuffle'; ...
                               'weighted Corr - unit ID shuffle' ; ...
                               'weighted Corr - wPBE time swap'};

sessionNames = {'Achilles-Maze1'; 'Achilles-Maze2'; 'Buddy'; 'Cicero'; 'Gatsby'; ...
                'Kevin';
                'RatN'; 'RatS'; 'RatU';
                'Roy-D1'; 'Roy-D2'; 'Roy-D3'; 'Ted-M1'; 'Ted-M2'; 'Ted-M3'};
                           
% initialize

nUnits = zeros(nSessions, 1);
peakPFlocation_session = cell(nSessions, 1);
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};

    % track spatial tunings
    spatialTuning_pooled = cell(nSessions, 1);
    
    
    % track learned tunings
    runLearnedTunings_pooled  = cell(nSessions, 1);
    
    
    % assembly tunings
    assemblyTuning_pooled.(currEpoch).data = cell(nSessions, 1);
    assemblyTuning_pooled.(currEpoch).ui   = cell(nSessions, 1);
    

    % PF correlation
    assemblyTuningPFcorr_pooled.(currEpoch).data   = cell(nSessions, 1);
    assemblyTuningPFcorr_pooled.(currEpoch).ui     = cell(nSessions, 1);

    assemblyTuningPFcorr_zscore_pooled.(currEpoch) = cell(nSessions, 1);

end

gw = gausswindow(5,15);

flag = 1;

% currSessions = 7:9;
currSessions = setdiff(1:15, 11:12);

for iSession = 1:numel(currSessions)
    
    currSession = currSessions(iSession);
    
    
    sessionName = sub(currSession+2).name
    
    
    load(fullfile(parentDir, sessionName, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr')
    spikes = spikes_pyr;
    
    s1 = load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'nPBEs', 'activeUnits', 'assemblyTunings', 'assemblyTuningPFcorr'); % , 'assemblyTunings', 'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore'
    for iEpoch = 1:3
        epochName = epochNames{iEpoch};
        assemblyTuningPFcorr.(epochName).ui = s1.assemblyTuningPFcorr.(epochName).ui;
        assemblyTuningPFcorr.(epochName).data = s1.assemblyTuningPFcorr.(epochName).data;
    end
    
    nPBEs = s1.nPBEs;
    activeUnits = s1.activeUnits;
    assemblyTunings = s1.assemblyTunings;
    
    
    
%     s2 = load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_1hour.mat']), 'assemblyTunings', 'assemblyTuningPFcorr');
%     for iEpoch = 1:3
%         epochName = epochNames{iEpoch};
%         assemblyTuningPFcorr.(epochName).data = s2.assemblyTuningPFcorr.(epochName).data;
%     end
%     assemblyTunings = s2.assemblyTunings;
    
    
    load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_activeRun.mat']), 'learnedTunings')
    
    
    okUnits = intersect(activeUnits.pre, activeUnits.post);
    okUnits = intersect(okUnits, activeUnits.run);
    
    spikes = spikes(okUnits);
    
    
    
    % track spatial tunings
    nUnits(iSession) = numel(spikes); 
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    
    
    interpPosBins = linspace(0, nPosBins, 200);
    nPosBins_interp = numel(interpPosBins);
    
    
    % non-directional spatial tuning
    spatialTunings_merge = zeros(nUnits(iSession), nPosBins_interp);
    
    peakPFlocation    = zeros(nUnits(iSession), 1);
    peakPFfiring      = zeros(nUnits(iSession), 1);
    
    
    
    for iUnit = 1:nUnits(iSession)
        
        currSpatialTuning = spikes(iUnit).spatialTuning_smoothed.uni;
        currSpatialTuning = interp1(1:nPosBins, currSpatialTuning, interpPosBins);
        currSpatialTuning = currSpatialTuning./max(currSpatialTuning);
        
        spatialTunings_merge(iUnit, :) = currSpatialTuning;
        [peakPFfiring(iUnit), peakPFlocation(iUnit)] = max(currSpatialTuning);
                
    end
    
    [~, sortIdx] = sort(peakPFlocation, 'ascend');
    sortedUnitIdx = okUnits(sortIdx);
    peakPFlocation_session{iSession} = peakPFlocation(sortIdx);
    
    
    temp = sum(spatialTunings_merge, 1);
    if flag == 1
        startBin = find(~isnan(temp), 1, 'first');
        endBin   = find(~isnan(temp), 1, 'last');
        flag = 0;
    end
    
    spatialTuning_pooled{iSession} = spatialTunings_merge(sortIdx, startBin:endBin);
    
    
    nPosBins_interp = endBin - startBin + 1;
    
        
    for iEpoch = 1:3
        currEpoch = epochNames{iEpoch};
        
        
        % assembly Tunings
        assemblyTuning_pooled.(currEpoch).data{iSession} = zeros(nUnits(iSession), nPosBins_interp);
        
        for iUnit = 1:nUnits(iSession)
            
            currLearnedTuning = assemblyTunings.(currEpoch).data(sortedUnitIdx(iUnit), :);
            currLearnedTuning = interp1(1:nPosBins, currLearnedTuning, interpPosBins);
            currLearnedTuning(isnan(currLearnedTuning)) = 0;
            currLearnedTuning = conv(currLearnedTuning, gw, 'same');
            currLearnedTuning = currLearnedTuning(startBin:endBin);
            
            assemblyTuning_pooled.(currEpoch).data{iSession}(iUnit, :) = currLearnedTuning/max(currLearnedTuning);
             
        end

        % PF correlation
        assemblyTuningPFcorr_pooled.(currEpoch).data{iSession}   = assemblyTuningPFcorr.(currEpoch).data(sortedUnitIdx);
        
        % PF correlation of learned tunings calcualted based on unit-identity shuffled PBEs 
        assemblyTuningPFcorr_pooled.(currEpoch).ui{iSession}   = assemblyTuningPFcorr.(currEpoch).ui(sortedUnitIdx, :);
    end
    
end

%%


cnctSpatialTunings = cell2mat(spatialTuning_pooled);  

cnctPFcorrelation = cell2mat([assemblyTuningPFcorr_pooled.pre.ui; assemblyTuningPFcorr_pooled.run.ui; assemblyTuningPFcorr_pooled.post.ui]);
% thresh = nanprctile(cnctPFcorrelation, 95);
thresh = -inf;


pre_hiIdx  = cell2mat(assemblyTuningPFcorr_pooled.pre.data) >= thresh;
run_hiIdx  = cell2mat(assemblyTuningPFcorr_pooled.run.data) >= thresh;
post_hiIdx = cell2mat(assemblyTuningPFcorr_pooled.post.data) >= thresh;



% z1   = numel(find(pre_hiIdx  & ~run_hiIdx & ~post_hiIdx));
% z2   = numel(find(~pre_hiIdx &  run_hiIdx & ~post_hiIdx));
% z3   = numel(find(~pre_hiIdx & ~run_hiIdx &  post_hiIdx));
% 
% z12  = numel(find(pre_hiIdx  &  run_hiIdx & ~post_hiIdx));
% z13  = numel(find(pre_hiIdx  & ~run_hiIdx &  post_hiIdx));
% z23  = numel(find(~pre_hiIdx &  run_hiIdx &  post_hiIdx));
% 
% z123 = numel(find(pre_hiIdx  &  run_hiIdx &  post_hiIdx));
% 
% 
% 
% 
% % venn diagram of units with significantly matching learned tunings during
% % each epoch
% 
% plotheight = 5;
% plotwidth  = 5;
% fontsize   = 6;
% 
% f= figure;
% clf(f);
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperSize', [plotwidth plotheight]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
% 
% 
% venn([z1 z2 z3 z12 z13 z23 z123], 'FaceColor',{'b','k', 'r'},'FaceAlpha',{0.5,0.5, 0.5},'EdgeColor','none')
% axis square
% legend('pre', 'run', 'post', 'location', 'eastout', 'box', 'off', 'fontsize', fontsize)
% axis off




hiLearnedTuning_unitIdx = find(pre_hiIdx | run_hiIdx | post_hiIdx);
hi_nUnits = numel(hiLearnedTuning_unitIdx);

loLearnedTuning_unitIdx = setdiff(1:sum(nUnits), hiLearnedTuning_unitIdx);
lo_nUnits = numel(loLearnedTuning_unitIdx);



plotheight = 11;
plotwidth  = 11;
fontsize   = 8;

nor = 1;
noc = 5;

leftmargin = 2;  rightmargin = 2;    topmargin = 2;    bottommargin = 2;

gapc = 0.2;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);




% for each panel

plotwidth_panel = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;

plotheight_hi = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr) * hi_nUnits/(hi_nUnits + lo_nUnits);
plotheight_lo = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr) * lo_nUnits/(hi_nUnits + lo_nUnits);


nor_panel = 2;
noc_panel = 1;

leftmargin_panel = 0;  rightmargin_panel = 0;     bottommargin_panel = 0;    topmargin_panel = 0;
gapr_panel = 0.5;    gapc_panel = 0.25;

% %

xfirst = leftmargin_panel;
yfirst = bottommargin_panel;
sub_pos_panel{1} = [xfirst yfirst plotwidth_panel plotheight_lo];


xfirst = leftmargin_panel;
yfirst = bottommargin_panel + plotheight_lo + gapr_panel;
sub_pos_panel{2} = [xfirst yfirst plotwidth_panel plotheight_hi];


sub_pos_panel_2 = sub_pos_panel;
for ii = 1 : 2
    sub_pos_panel_2{ii}(1) = (sub_pos_panel{ii}(1) + sub_pos{1}(1)*plotwidth) / plotwidth;
    sub_pos_panel_2{ii}(2) = (sub_pos_panel{ii}(2) + sub_pos{1}(2)*plotheight) / plotheight;
    sub_pos_panel_2{ii}(3) = sub_pos_panel{ii}(3) / plotwidth;
    sub_pos_panel_2{ii}(4) = sub_pos_panel{ii}(4) / plotheight;
end




% % place fields

% best learned tunings

hiSpatialTunings = cnctSpatialTunings(hiLearnedTuning_unitIdx, :);

[~, peakLocs] = max(hiSpatialTunings, [], 2);
[~, hiSortIdx]  = sort(peakLocs, 'ascend');

hiSpatialTunings = hiSpatialTunings(hiSortIdx, :);



% worst learned tunings

loSpatialTunings = cnctSpatialTunings(loLearnedTuning_unitIdx, :);

[~, peakLocs] = max(loSpatialTunings, [], 2);
[~, loSortIdx]  = sort(peakLocs, 'ascend');

loSpatialTunings = loSpatialTunings(loSortIdx, :);


% ax1 = axes('position',sub_pos_panel_2{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
% 
% imagesc(loSpatialTunings)
% xlim([0 nPosBins_interp])
% ylim([1 lo_nUnits])
% set(ax1, 'xTick', [], 'yTick', [1 lo_nUnits], 'YDir', 'normal', 'box', 'off')
% 
% h = text(-200, 0, 'LT-PF corr < 95% of null', 'fontsize', fontsize);
% set(h, 'rotation', 90)
% 
% h = text(-100, lo_nUnits/2, 'units', 'fontsize', fontsize);
% set(h, 'rotation', 90)




ax2 = axes('position',sub_pos_panel_2{2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

imagesc(hiSpatialTunings)
xlim([0 nPosBins_interp])
ylim([1 hi_nUnits])
set(ax2, 'xTick', [0 nPosBins_interp], 'xTickLabel', {'0'; '1'}, 'yTick', [1 hi_nUnits], 'YDir', 'normal', 'box', 'off')

h = text(-200, 0, 'LT-PF corr > 95% of null', 'fontsize', fontsize);
set(h, 'rotation', 90)

h = text(-100, hi_nUnits/2, 'units', 'fontsize', fontsize);
set(h, 'rotation', 90)

dim = [sub_pos_panel_2{2}(1)-0.02 sub_pos_panel_2{2}(2)+0.12 0.3 0.3];
annotation('textbox', dim,'String','place fields' ,'FitBoxToText','on', 'EdgeColor', 'none', 'fontsize', fontsize)
 



for iEpoch = 1:3

    epochName = epochNames{iEpoch};
        
    cnctAssemblyTunings = cell2mat(assemblyTuning_pooled.(epochName).data);
    

    % best learned tunings
    hiLearnedTunings = cnctAssemblyTunings(hiLearnedTuning_unitIdx, :);
    hiLearnedTunings = hiLearnedTunings(hiSortIdx, :);
    

    % worst learned tunings
    loLearnedTunings = cnctAssemblyTunings(loLearnedTuning_unitIdx, :);    
    loLearnedTunings = loLearnedTunings(loSortIdx, :);

    
    
    sub_pos_panel_2 = sub_pos_panel;
    for ii = 1 : 2
        sub_pos_panel_2{ii}(1) = (sub_pos_panel{ii}(1) + sub_pos{iEpoch+2}(1)*plotwidth) / plotwidth;
        sub_pos_panel_2{ii}(2) = (sub_pos_panel{ii}(2) + sub_pos{iEpoch+2}(2)*plotheight) / plotheight;
        sub_pos_panel_2{ii}(3) = sub_pos_panel{ii}(3) / plotwidth;
        sub_pos_panel_2{ii}(4) = sub_pos_panel{ii}(4) / plotheight;
    end
    
%     ax1 = axes('position',sub_pos_panel_2{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
%     imagesc(loLearnedTunings)
%     xlim([0 nPosBins_interp])
%     ylim([1 lo_nUnits])
%         
%     set(ax1, 'xTick', [0 nPosBins_interp], 'xTickLabel', {'0'; '1'}, 'yTick', [], 'YDir', 'normal', 'box', 'off')

    

    
    if iEpoch == 2
        xlabel('position(normalized)', 'fontsize', fontsize)
    end
    
    
    ax2 = axes('position',sub_pos_panel_2{2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    imagesc(hiLearnedTunings)
    xlim([0 nPosBins_interp])
    ylim([1 hi_nUnits])
    
    set(ax2, 'xTick', [], 'yTick', [], 'YDir', 'normal', 'box', 'off')


    
    dim = [sub_pos_panel_2{2}(1)-0.02 sub_pos_panel_2{2}(2)+0.12 0.3 0.3];
    annotation('textbox',dim,'String',[epochName ' PBEs'] ,'FitBoxToText','on', 'EdgeColor', 'none', 'fontsize', fontsize)
    
    
    
    
end
dim = [sub_pos{3}(1)-0.05 sub_pos_panel_2{2}(2)+0.2 0.3 0.3];
annotation('textbox',dim,'String','learned tunings' ,'FitBoxToText','on', 'EdgeColor', 'none')
    


colormap('jet')


%% functions

function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end

