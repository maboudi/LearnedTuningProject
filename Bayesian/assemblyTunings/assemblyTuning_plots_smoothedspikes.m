clear
clc


sessionNumbers = [14:15];

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults';

files = dir(parentDir);
isdir = [files.isdir];
subfolders = files(isdir);

sessionNames = {subfolders(3:end).name};

epochNames = {'pre'; 'run'; 'post'};

faceColors = {'b'; 'none';'r'};
edgeColors = {'none'; 'k'; 'none'};
colors     = {'b'; 'k';'r'};

Qs = (0:5) * 0.02; % in sec


%% preconfigureing the plots


plotheight = 27.5; %% in cm
plotwidth = 20;
fontsize = 6;

nor = 3;
noc = 2; %% number of events per column

leftmargin = 2;  rightmargin = 2;    topmargin = 4;    bottommargin = 2;

gapc = 3;
gapr = 2.5;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f=figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
    



plotwidth2 = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
plotheight2 = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr)/nor;

nor2 = 2;
noc2 = 3;

leftmargin2 = 0;  rightmargin2 = 0;     bottommargin2 = 0;    topmargin2 = 0;
gapr2 = 1.5;    gapc2 = 1;

sub_pos2 = subplot_pos(plotwidth2, plotheight2, leftmargin2, rightmargin2, bottommargin2, topmargin2, noc2, nor2, gapr2, gapc2); %% this gives us the positions for the subplots, each corresponded to an event



for isession = 1:numel(sessionNumbers)
    
   sessionNumber = sessionNumbers(isession); 
   sessionName   = sessionNames{sessionNumber};
   
   basePath = fullfile(parentDir, sessionName);

   
   
   load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']))
   asTuning_PF_spatailBinCorr_20ms    = assemblyTuningPFcorr;
   asTuning_PF_prefPosDist_20ms       = asTuning_PF_prefPosDist;
   asTuning_PF_prefPos_spearCorr_20ms = asTuning_PF_prefPos_spearCorr;
   asTuning_PF_prefPos_pval_20ms      = asTuning_PF_prefPos_pval;
   
   
   load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_50ms.mat']))
   asTuning_PF_spatailBinCorr_50ms    = assemblyTuningPFcorr;
   asTuning_PF_prefPosDist_50ms       = asTuning_PF_prefPosDist;
   asTuning_PF_prefPos_spearCorr_50ms = asTuning_PF_prefPos_spearCorr;
   asTuning_PF_prefPos_pval_50ms      = asTuning_PF_prefPos_pval;
   
   
   load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_smoothedPBEs.mat']))
   asTuning_PF_spatailBinCorr_smoothed    = assemblyTuningPFcorr;
   asTuning_PF_prefPosDist_smoothed       = asTuning_PF_prefPosDist;
   asTuning_PF_prefPos_spearCorr_smoothed = asTuning_PF_prefPos_spearCorr;
   asTuning_PF_prefPos_pval_smoothed      = asTuning_PF_prefPos_pval;
   
   
   for iepoch = 1:3
       
        currEpoch = epochNames{iepoch};
       
        sub_pos2 = subplot_pos(plotwidth2, plotheight2, leftmargin2, rightmargin2, bottommargin2, topmargin2, noc2, nor2, gapr2, gapc2); 
        for ii = 1 : nor2
            for jj = 1 : noc2
                sub_pos2{ii,jj}(1) = (sub_pos2{ii,jj}(1) * plotwidth2 + sub_pos{iepoch, isession}(1) * plotwidth) / plotwidth;
                sub_pos2{ii,jj}(2) = (sub_pos2{ii,jj}(2) * plotheight2 + sub_pos{iepoch, isession}(2) * plotheight) / plotheight;
                sub_pos2{ii,jj}(3) = sub_pos2{ii,jj}(3) * plotwidth2 / plotwidth;
                sub_pos2{ii,jj}(4) = sub_pos2{ii,jj}(4) * plotheight2 / plotheight;
            end
        end
        
        
        % plot spatial bin correlation for different smoothing Qs
        axes('position',sub_pos2{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');

        h = violinplot(asTuning_PF_spatailBinCorr_smoothed.(currEpoch));

        for ii = 1:numel(Qs)
            h(ii).ViolinColor = faceColors{iepoch};
            h(ii).EdgeColor   = edgeColors{iepoch};
            h(ii).ShowData    = 0;
            h(ii).MedianPlot.SizeData = 10;
            h(ii).MedianPlot.Marker = 'square';
            h(ii).MedianPlot.MarkerEdgeColor = 'k';
            h(ii).MedianPlot.MarkerFaceColor = 'none';
            h(ii).WhiskerPlot.LineWidth = 1;
        end
        
        xticklabels({'0'; '.02'; '.04'; '.06'; '.08'; '.1'})
        xtickangle(45)
        xlabel('smoothing SD', 'fontsize', fontsize) 
        
        ylabel('AT-PF spatial bin correlation', 'fontsize', fontsize)
        ylim([-1 1])
        grid on
        
        
        
        
        % plot preferred location distance for different smoothing Qs
        axes('position',sub_pos2{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');

        h = violinplot(asTuning_PF_prefPosDist_smoothed.(currEpoch));

        for ii = 1:numel(Qs)
            h(ii).ViolinColor = faceColors{iepoch};
            h(ii).EdgeColor   = edgeColors{iepoch};
            h(ii).ShowData    = 0;
            h(ii).MedianPlot.SizeData = 10;
            h(ii).MedianPlot.Marker = 'square';
            h(ii).MedianPlot.MarkerEdgeColor = 'k';
            h(ii).MedianPlot.MarkerFaceColor = 'none';
            h(ii).WhiskerPlot.LineWidth = 1;
        end
        
        xticklabels({'0'; '.02'; '.04'; '.06'; '.08'; '.1'})
        xtickangle(45)
        xlabel('smoothing SD', 'fontsize', fontsize)
        
        ylabel('AT-PF pref. location distance', 'fontsize', fontsize)
        ylim([0 1])
        grid on
        
        
        % plot preferred location correlation for different smoothing Qs
        axes('position',sub_pos2{1,3},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');
        
        plot(1:numel(Qs), asTuning_PF_prefPos_spearCorr_smoothed.(currEpoch), 'marker', 'square', 'MarkerSize', 3, 'color', colors{iepoch});
        
        for iq = 1:numel(Qs)
            text(iq, 0.8+0.025, calSignificance_sign(asTuning_PF_prefPos_pval_smoothed.(currEpoch)(iq)), 'fontsize', fontsize, 'color', 'k')
        end

        xticks(1:numel(Qs))
        xticklabels({'0'; '.02'; '.04'; '.06'; '.08'; '.1'})
        xtickangle(45)
        xlabel('smoothing SD', 'fontsize', fontsize) 
        
        
        ylabel('AT-PF pref. location correlation', 'fontsize', fontsize)
        ylim([-0.8 0.8])
        grid on
        
        set(gca, 'box', 'off', 'FontSize',fontsize)
        
        
        
        
        
        % plot AT-PF spatial bin correlation 20 vs 50 ms time bins
        axes('position',sub_pos2{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');

        h = violinplot([asTuning_PF_spatailBinCorr_20ms.(currEpoch).data asTuning_PF_spatailBinCorr_50ms.(currEpoch).data]);

        for ii = 1:2
            h(ii).ViolinColor = faceColors{iepoch};
            h(ii).EdgeColor   = edgeColors{iepoch};
            h(ii).ShowData    = 0;
            h(ii).MedianPlot.SizeData = 10;
            h(ii).MedianPlot.Marker = 'square';
            h(ii).MedianPlot.MarkerEdgeColor = 'k';
            h(ii).MedianPlot.MarkerFaceColor = 'none';
            h(ii).WhiskerPlot.LineWidth = 1;
        end
        
        xlim([-1 4])
        xticklabels({'.02'; '.05'})
        xtickangle(45)
        xlabel('time bin size', 'fontsize', fontsize) 
        
        ylabel('AT-PF spatial bin correlation', 'fontsize', fontsize)
        ylim([-1 1])
        grid on
        
        
        
        
        
        % plot preferred location distance for different smoothing Qs
        axes('position',sub_pos2{2,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');

        h = violinplot([asTuning_PF_prefPosDist_20ms.(currEpoch).data asTuning_PF_prefPosDist_50ms.(currEpoch).data]);

        for ii = 1:2
            h(ii).ViolinColor = faceColors{iepoch};
            h(ii).EdgeColor   = edgeColors{iepoch};
            h(ii).ShowData    = 0;
            h(ii).MedianPlot.SizeData = 10;
            h(ii).MedianPlot.Marker = 'square';
            h(ii).MedianPlot.MarkerEdgeColor = 'k';
            h(ii).MedianPlot.MarkerFaceColor = 'none';
            h(ii).WhiskerPlot.LineWidth = 1;
        end
        
        xlim([-1 4])
        xticklabels({'.02'; '.05'})
        xtickangle(45)
        xlabel('time bin size', 'fontsize', fontsize) 

        ylabel('AT-PF pref. location distance', 'fontsize', fontsize)
        ylim([0 1])
        grid on 
        
        
        % plot preferred location correlation for different smoothing Qs
        axes('position',sub_pos2{2,3},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');
        
        
        currData = [];
        
        plot(1:2, [asTuning_PF_prefPos_spearCorr_20ms.(currEpoch).data asTuning_PF_prefPos_spearCorr_50ms.(currEpoch).data], 'marker', 'square', 'MarkerSize', 3, 'color', colors{iepoch});
        
        text(1, 0.8+0.025, calSignificance_sign(asTuning_PF_prefPos_pval_20ms.(currEpoch).data), 'fontsize', fontsize, 'color', 'k')
        text(2, 0.8+0.025, calSignificance_sign(asTuning_PF_prefPos_pval_50ms.(currEpoch).data), 'fontsize', fontsize, 'color', 'k')
        
        xlim([-1 4])
        xticks(1:2)
        xticklabels({'.02'; '.05'}) 
        xtickangle(45)
        xlabel('time bin size', 'fontsize', fontsize)
        
        
        ylabel('AT-PF pref. location correlation', 'fontsize', fontsize)
        ylim([-0.8 0.8])
        grid on
        
        set(gca, 'box', 'off', 'FontSize',fontsize)
        
   end

end

sessionName = sessionNames{sessionNumbers(1)};
sessionName(sessionName == '_') = '-';
annotation('textbox', [0.075 0.85 0.1 0.1], 'string', sessionName, 'fontsize', 10, 'EdgeColor', 'none')

annotation('textbox', [0.075 0.80 0.1 0.1], 'string', 'PRE', 'fontsize', 10, 'EdgeColor', 'none', 'fontweight', 'bold')
annotation('textbox', [0.075 0.50 0.1 0.1], 'string', 'RUN', 'fontsize', 10, 'EdgeColor', 'none', 'fontweight', 'bold')
annotation('textbox', [0.075 0.21 0.1 0.1], 'string', 'POST', 'fontsize', 10, 'EdgeColor', 'none', 'fontweight', 'bold')


if numel(sessionNumbers) == 2
    sessionName = sessionNames{sessionNumbers(2)};
    sessionName(sessionName == '_') = '-';
    annotation('textbox', [0.55 0.85 0.1 0.1], 'string', sessionName, 'fontsize', 10, 'EdgeColor', 'none')

    annotation('textbox', [0.55 0.80 0.1 0.1], 'string', 'PRE', 'fontsize', 10, 'EdgeColor', 'none', 'fontweight', 'bold')
    annotation('textbox', [0.55 0.50 0.1 0.1], 'string', 'RUN', 'fontsize', 10, 'EdgeColor', 'none', 'fontweight', 'bold')
    annotation('textbox', [0.55 0.21 0.1 0.1], 'string', 'POST', 'fontsize', 10, 'EdgeColor', 'none', 'fontweight', 'bold')
end

print(gcf, 'set_8.pdf', '-dpdf', '-painters')
close all


%% functions


function significance_sign = calSignificance_sign(pvalue)    

if pvalue < 0.001
    significance_sign = '***';
elseif pvalue < 0.01
    significance_sign = '**';
elseif pvalue < 0.05
    significance_sign = '*';
else
    significance_sign = '';
end


end
