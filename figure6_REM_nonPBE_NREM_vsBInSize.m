
clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 10 13:15 6];
nSessions    = numel(currSessions);

binDur = [0.02 0.05 0.125 0.25];

epochNames = {'pre'; 'run'; 'post'};
dataTypes  = {'data'; 'ui'};


for iEpoch = [1 3]
    currEpoch = epochNames{iEpoch};
    
    for idata = 1:2
        currDataType = dataTypes{idata};
        learnedTuningPFcorr.(currEpoch).(currDataType)  = cell(nSessions, 1);
        learnedTuningPFKLdiv.(currEpoch).(currDataType) = cell(nSessions, 1);
    end
    
    PFcorr_ranksum_pval.(currEpoch)  = nan(3,1);
    PFKLdiv_ranksum_pval.(currEpoch) = nan(3,1);
    
end



tt = 1;

for iBin = 1:numel(binDur)
    
    currBinDur = binDur(iBin);

    for iSess = 1:nSessions

        sessionNumber = currSessions(iSess);
        sessionName   = rr(sessionNumber+2).name;
        basePath      = fullfile(parentDir, sessionName);


        % % maybe for now going ahead with the active units that were calculated
        % for the PBEs

        s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), 'activeUnits');

        activeUnits = s.activeUnits;
        okUnits = intersect(activeUnits.pre, activeUnits.post);
        okUnits = intersect(okUnits, activeUnits.run);


        % % learned tunings and PF-matching scores
        s = load(fullfile(basePath, 'assemblyTunings', 'NREM_REM', sprintf('%s.assemblyTunings_REM_%.3f.mat', sessionName, currBinDur)), ...
            'assemblyTuningPFcorr', 'assemblyTuningPFKLdiv');


        for iEpoch = [1 3]
            currEpoch = epochNames{iEpoch};

            for idata = 1:2
                currDataType = dataTypes{idata};

                learnedTuningPFcorr.(currEpoch).(currDataType){iSess}  = s.assemblyTuningPFcorr.(currEpoch).(currDataType)(okUnits, :);
                learnedTuningPFKLdiv.(currEpoch).(currDataType){iSess}  = s.assemblyTuningPFKLdiv.(currEpoch).(currDataType)(okUnits, :);
            end

        end

    end


    for iEpoch = [1 3]
        currEpoch = epochNames{iEpoch};
        
        currPFcorrs = learnedTuningPFcorr.(currEpoch);
        currPFKLdiv = learnedTuningPFKLdiv.(currEpoch);
        
        
        for idata = 1:2
            currDataType = dataTypes{idata};

            currPFcorrs.(currDataType)  = cell2mat(currPFcorrs.(currDataType));
            currPFKLdiv.(currDataType)  = cell2mat(currPFKLdiv.(currDataType));

            PFcorr_stats.(currEpoch).(currDataType)(tt, :)  = [nanprctile(currPFcorrs.(currDataType), 25) nanmedian(currPFcorrs.(currDataType)) nanprctile(currPFcorrs.(currDataType), 75)];
            PFKLdiv_stats.(currEpoch).(currDataType)(tt, :) = [nanprctile(currPFKLdiv.(currDataType), 25) nanmedian(currPFKLdiv.(currDataType)) nanprctile(currPFKLdiv.(currDataType), 75)];  
        end
        
%         PFcorr_ranksum_pval.(currEpoch)(tt)  = ranksum(currPFcorrs.data, currPFcorrs.ui, 'tail', 'right');
%         PFKLdiv_ranksum_pval.(currEpoch)(tt) = ranksum(currPFKLdiv.data, currPFKLdiv.ui, 'tail', 'left');
        
%         PFcorr_ranksum_pval.(currEpoch)(tt)  = signrank(currPFcorrs.data, 0, 'tail', 'right');
%         PFKLdiv_ranksum_pval.(currEpoch)(tt) = ranksum(currPFKLdiv.data, currPFKLdiv.ui, 'tail', 'left');
        
 
    end

    tt = tt + 1;

end


%%
% plot the distribution of LT-PFcorrelation


currVariable = PFKLdiv_stats;

plotheight = 61;
plotwidth  = 60;
fontsize   =  4;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

colors = [[0 0 200]/255; [0 153 0]/255; [204 0 0]/255];


for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
    currStats = currVariable.(currEpoch).data;
    plot((1:numel(binDur))+0.07*iEpoch, currStats(:,2), 'color', [colors(iEpoch, :) 0.5], 'linewidth', 1)

    for iBin = 1:numel(binDur)
        
        currBinDur = binDur(iBin);
        line([iBin iBin]+0.07*iEpoch, currStats(iBin, [1 3]), 'color', [colors(iEpoch, :) 0.5], 'linewidth', 1)

    end
    
end


xlim([0.5 4.5])
% ylim([-0.4 0.4])
yl = ylim;
ylim([0 3])

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])
set(gca, 'XTick', 1:4, 'XTickLabel', {'0.02'; '0.05'; '0.125'; '0.25'})
% set(gca, 'yTick', -0.4:0.2:0.4)
set(gca, 'yTick', 0:0.5:3)

xtickangle(gca, 45)
aa = gca;
aa.YGrid = 'on';

ylabel('KL divergence', 'fontsize', fontsize)

% ylabel('LT-PF KL div.', 'fontsize', fontsize)
xlabel('time bin (sec)', 'fontsize', fontsize)

% axis square



%% sub-functions

function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end

