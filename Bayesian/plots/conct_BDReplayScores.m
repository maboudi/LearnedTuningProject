clear; clc; close all

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

sub = dir(parentDir);
% nSessions = numel(sub) - 2; % subtracting two becuase of '.' and '..'


epochNames = {'pre'; 'run'; 'post'};

replayScoreMethods = {'rt_ui'; 'wc_ts'; 'wc_pf'; 'wc_ui'; 'wc_ds'};


replayScoreMethods_fullName = {'radon Integral - unit ID shuffle'; ...
                               'weighted Corr - wPBE time swap'};


currSessions = [1:5 7:9 12 15:17 6];
nSessions    = numel(currSessions);
                           
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};
    
    for irsm = 1:numel(replayScoreMethods)
        replayScores.(replayScoreMethods{irsm}).(currEpoch) = cell(nSessions, 1);   
    end
end

for iSess = 1:nSessions


    currSession = currSessions(iSess);

    sessionName = sub(currSession+2).name


    load(fullfile(parentDir, sessionName, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
    behavior = fileInfo.behavior;

    startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
    startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
    startT.post = behavior.time(2,2); endT.post = startT.post + 4*60*60; %behavior.time(3,2); 


    
    load(fullfile(parentDir, sessionName, 'BayesianDecoding', [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')
    
    
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};
        
        idx = [PBEInfo_replayScores.peakT] > startT.(currEpoch) & [PBEInfo_replayScores.peakT] < endT.(currEpoch);


        for irsm = 1:numel(replayScoreMethods)

            replayScores.(replayScoreMethods{irsm}).(currEpoch){iSess} = ...
                [PBEInfo_replayScores(idx).(replayScoreMethods{irsm})];

        end

        replayScores.wc.(currEpoch){iSess} = [PBEInfo_replayScores(idx).weightedCorr];
        replayScores.radon.(currEpoch){iSess} = [PBEInfo_replayScores(idx).radonIntegral];

    end

end


%% concatenating the results from differnet sessions to plot the distributions


currSessions = 6:8;

for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};

    for irsm = 1:numel(replayScoreMethods)

        replayScores_pooled.(replayScoreMethods{irsm}).(currEpoch) = ...
            cell2mat(replayScores.(replayScoreMethods{irsm}).(currEpoch)(currSessions)');
    end

    replayScores_pooled.wc.(currEpoch) = cell2mat(replayScores.wc.(currEpoch)(currSessions));
    replayScores_pooled.radon.(currEpoch) = cell2mat(replayScores.radon.(currEpoch)(currSessions));

end


%%

plotwidth  = 150;
plotheight = 80;
fontsize   = 6;

f = figure;

clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


% spatial bin correlation
customCDF(replayScores_pooled.wc)

xlabel({'replay score %'}, 'fontsize', fontsize)
ylabel('fraction of ripple events', 'fontsize', fontsize)

% rr = legend(gca, 'pre', 'run', 'post', 'box', 'off', 'location', 'northout');
% rr.FontSize = fontsize;

set(gca, 'box', 'off','linewidth', 1, 'TickDir', 'out')



%% subfunctions

function customCDF(variable)


epochNames = {'pre'; 'run'; 'post'};

if isfield(variable.pre, 'data') % in case we are plotting the absolute values and not z-scores
    
    for iEpoch = 1:3
        variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
    end
    
    variable = variable2;
    clear variabale2
end 


colors = [[0 92 233]/255; [0 0 0]; [212 0 53]/255];

hold on
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};
    pooledVar = variable.(currEpoch);     
    
    nUnits = numel(pooledVar);

    [cdf_pooled, x_pooled] = ecdf(pooledVar);
    auc.(currEpoch) = trapz(x_pooled, cdf_pooled);

    [x_pooled, idx] = unique(x_pooled);
    cdf_pooled = cdf_pooled(idx);
    
%     area(x_pooled, cdf_pooled, 'FaceColor', colors(iEpoch, :), 'FaceAlpha', 0.3, 'EdgeColor', colors(iEpoch, :), 'linewidth', 0.5, 'EdgeAlpha', 0.6)
    plot(x_pooled, cdf_pooled, 'color', colors(iEpoch, :), 'linewidth', 1)

    medianCorr.(currEpoch) = nanmedian(pooledVar);

end


nIters = 10000;
auc_chance = nan(nIters, 1);


for iIter = 1:nIters
    rndPrctls = randi(100, nUnits, 1);

    [cdf_chnace, x_chance] = ecdf(rndPrctls);

    auc_chance(iIter) = trapz(x_chance, cdf_chnace);
end



ylim([0 1])
yl = ylim;
xl = xlim;


for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    [~, pval] = ttest2(auc.(currEpoch), auc_chance, 'tail', 'left');

    
    line([xl(1) medianCorr.(currEpoch)], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors(iEpoch, :))
    line([medianCorr.(currEpoch) medianCorr.(currEpoch)], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors(iEpoch, :))
    text(medianCorr.(currEpoch), 0.07, sprintf('%.2f(p=%.2e)', medianCorr.(currEpoch), pval), 'fontsize', 6, 'color', colors(iEpoch, :))
end


% significance scores
ranksumPvalue = nan(3);
for iEpoch = 1:3
   var1 = variable.(epochNames{iEpoch});
   
   for jEpoch = setdiff(1:3, iEpoch)
       var2 = variable.(epochNames{jEpoch});

       ranksumPvalue(iEpoch, jEpoch) = ranksum(var1, var2, 'tail', 'right');
   end
end
       
% add the significance signs here ...  

xticks([0:25:100])
yticks([0:0.5:1])
    
end








