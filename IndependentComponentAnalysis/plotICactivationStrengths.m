function plotICactivationStrengths(activationStrength, binCenters, behavior, nor, noc, timeScale)


periods      = {'PRE'; 'RUN'; 'POST'};
% periodColors = [[51 150 255]/255; [44 47 50]/255; [220 90 90]/255];
periodColors = [[220 90 90]/255; [220 90 90]/255; [220 90 90]/255];


states      = {'NREM'; 'Wake'; 'REM'; 'Quiet'};
stateColors = [[0 0 1]; [1 1 1]; [1 0 0]; [0 0 0]];

state = struct('NREM', [], 'Wake', [], 'REM', [], 'Quiet', []);
state.NREM  = changeFormat([behavior.NREM; behavior.Intermediate], behavior);
state.Wake  = changeFormat(behavior.Wake, behavior);
state.REM   = changeFormat(behavior.REM, behavior);
state.Quiet = changeFormat(behavior.Drowsy, behavior);

binCenters.PRE  = binCenters.PRE  - behavior.time(1,1);
binCenters.RUN  = binCenters.RUN  - behavior.time(2,1);
binCenters.POST = binCenters.POST - behavior.time(3,1);


%%% configure the plots

plotheight = 40; %% in cm
plotwidth  = 50;

% nor = 1; %% number of events per row
% noc = 1; %% number of events per column

% leftedge = 10;
% rightedge = 10;
% topedge = 10;
% bottomedge = 10;


leftedge = 2;
rightedge = 2;
topedge = 5;
bottomedge = 2;


gapr = 2;
gapc = 3;

fontsize = 7;


%%% each subplot (corresponding to an pattern) is comprised of 3 panels corresponding to PRE, RUN, and POST reactivation strength time courses 

plotwidth_inset  = (plotwidth - (rightedge+leftedge) - (noc-1) * gapc)/noc;
plotheight_inset = (plotheight - (topedge+bottomedge) - (nor-1) * gapr)/nor;

nor_inset = 1;
noc_inset = 3;

leftedge_inset = 0;
rightedge_inset = 0;
bottomedge_inset = 0;
topedge_inset = 0;

gapr_inset = 0.2;
gapc_inset = 0;


sub_pos  = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); %% this gives us the positions for the subplots, each corresponded to an event


f=figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

flag = 0;
    
for i = 1:nor
    
    if flag == 1
        break
    end
    
    for j = 1:noc
        patternInd = (i - 1) * noc + j;
        
        if patternInd > size(activationStrength.PRE.original, 1)
            flag = flag+1;
            break
        end
        

        sub_pos_inset = subplot_pos_inrow2(plotwidth_inset, plotheight_inset, leftedge_inset, rightedge_inset, bottomedge_inset, topedge_inset, noc_inset, nor_inset, gapr_inset, gapc_inset, [1 numel(binCenters.RUN)/numel(binCenters.PRE) numel(binCenters.POST)/numel(binCenters.PRE)]); %% this is a modified version of subplot_pos


        for ii = 1 : size(sub_pos_inset, 1)
            for jj = 1 : noc_inset

                sub_pos_inset{ii,jj}(1) = (sub_pos_inset{ii,jj}(1) * plotwidth_inset + sub_pos{i,j}(1) * plotwidth) / plotwidth;
                sub_pos_inset{ii,jj}(2) = (sub_pos_inset{ii,jj}(2) * plotheight_inset + sub_pos{i,j}(2) * plotheight) / plotheight;
                sub_pos_inset{ii,jj}(3) = sub_pos_inset{ii,jj}(3) * plotwidth_inset / plotwidth;
                sub_pos_inset{ii,jj}(4) = sub_pos_inset{ii,jj}(4) * plotheight_inset / plotheight;
            end
        end
        
        
        pooledData = [activationStrength.PRE.original(patternInd, :) ...
            activationStrength.POST.original(patternInd, :)];
%                       activationStrength.RUN.original(patternInd, :) ...
                      
        
        pooledDataSmoothed = [activationStrength.PRE.smoothed(patternInd, :) ...
                              activationStrength.POST.smoothed(patternInd, :)];
                  
        
        for prdIdx = [1 3] % period index (PRE, RUN, POST)
            
            axes('position',sub_pos_inset{1, prdIdx},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
            hold on
            set(gca, 'linewidth',1.5, 'TickDir', 'out')
            
            hAx = plotyy(binCenters.(periods{prdIdx})/3600, activationStrength.(periods{prdIdx}).original(patternInd, :), ...
                binCenters.(periods{prdIdx})/3600, activationStrength.(periods{prdIdx}).smoothed(patternInd, :));
            
%             hAx = plotyy(binCenters.(periods{prdIdx}), activationStrength.(periods{prdIdx}).original(patternInd, :), ...
%                 binCenters.(periods{prdIdx}), activationStrength.(periods{prdIdx}).smoothed(patternInd, :));
            
            
            hAx(1).Children.Color     = [0.5 0.5 0.5];
            hAx(1).Children.LineWidth = 1; 
            
            if max(pooledData > 200)
                hAx(1).YTick = 0:100:max(pooledData);
            else
                hAx(1).YTick = 0:50:max(pooledData);
            end
            
%             hAx(1).YLim               = [min(pooledData) max(pooledData)];
            hAx(1).YAxis.Color        = [0.15 0.15 0.15];
            hAx(1).XLim               = (behavior.time(prdIdx, :) - behavior.time(prdIdx, 1))/3600;
            
            hAx(2).Children.Color     = periodColors(prdIdx, :);
            hAx(2).Children.LineWidth = 2; 
            hAx(2).YTick              = [-2:2:max(2, max(pooledDataSmoothed))];
%             hAx(2).YLim               = [min(pooledDataSmoothed) max(pooledDataSmoothed)];
            hAx(2).YAxis.Color        = periodColors(prdIdx, :);
            hAx(2).XLim               = hAx(1).XLim;
            
            
            % if the diagonal of projection operator was set to zero(hence,
            % negative reactivations) and we are using z-scored
            % reactivation strength
            
            if abs(max(pooledData)/min(pooledData)) > abs(max(pooledDataSmoothed)/min(pooledDataSmoothed))
               hAx(1).YLim = [min(pooledData) max(pooledData)];
               hAx(2).YLim = [min(pooledDataSmoothed) min(pooledDataSmoothed)* max(pooledData)/min(pooledData)];
            else
               hAx(1).YLim = [min(pooledData) min(pooledData)*max(pooledDataSmoothed)/min(pooledDataSmoothed)]; 
               hAx(2).YLim = [min(pooledDataSmoothed) max(pooledDataSmoothed)];
            end
              
              % if the diagonal of projection operator was remained intact and not using z-scored reactivation strength 
%             hAx(1).YLim = [min(pooledData) max(pooledData)]; 
%             hAx(2).YLim = [min(pooledDataSmoothed) max(pooledDataSmoothed)];
            
            % if the diagonal of projection operator was remianed intact
            % and using z-scred reactivation strength
            
%             hAx(1).YLim = [max(pooledData)*min(pooledDataSmoothed)/max(pooledDataSmoothed) max(pooledData)]; 
%             hAx(2).YLim = [min(pooledDataSmoothed) max(pooledDataSmoothed)];
            
            
            dummy = ylim;
            for stateIdx = 1:4 % NREM, Wake, REM, Quiet
                
                currSegments = state.(states{stateIdx}).(periods{prdIdx}); % segments are each single RUN, NREM, Wake, or Quiet bouts/episodes
                
                if ~isempty(currSegments)
                    for iseg = 1: size(currSegments, 1)
                       patch([currSegments(iseg, 1) currSegments(iseg, 2) currSegments(iseg, 2) currSegments(iseg, 1)]/3600, [dummy(1) dummy(1) dummy(2) dummy(2)], stateColors(stateIdx, :), 'facealpha', 0.2, 'edgeColor', 'none')
%                        patch([currSegments(iseg, 1) currSegments(iseg, 2) currSegments(iseg, 2) currSegments(iseg, 1)], [dummy(1) dummy(1) dummy(2) dummy(2)], stateColors(stateIdx, :), 'facealpha', 0.2, 'edgeColor', 'none')

                    end
                end
            end
            
            xlim((behavior.time(prdIdx, :) - behavior.time(prdIdx, 1))/3600)
%             xlim((behavior.time(prdIdx, :) - behavior.time(prdIdx, 1)))

            xlabel('time(hr)', 'fontsize', fontsize)
            
            set(gca, 'box', 'off', 'xtick', [0:1:max(binCenters.(periods{prdIdx})/3600)])
%             set(gca, 'box', 'off', 'xtick', [0:1:max(binCenters.(periods{prdIdx}))])


            if prdIdx == 1
               hAx(1).YLabel.String   = 'reactivation strength';
               hAx(1).YLabel.FontSize = fontsize;
               hAx(2).YTickLabel      = [];
               title('PRE', 'fontsize', fontsize, 'fontweight', 'normal')
               text(0, 1.1*range(dummy)+dummy(1), sprintf('PC %d', patternInd), 'fontsize', 10)
            end
            
           if prdIdx == 2
              hAx(1).YTickLabel  = [];
              hAx(2).YTickLabel  = [];
              title('RUN', 'fontsize', fontsize, 'fontweight', 'normal');
            end
            
            if prdIdx == 3  
               hAx(1).YTickLabel      = [];
               hAx(2).YLabel.String   = 'average reactivation';
               hAx(2).YLabel.FontSize = fontsize;
               hAx(2).YLabel.Color    = periodColors(prdIdx, :);
               title('POST', 'fontsize', fontsize, 'fontweight', 'normal'); 
            end
        end
            
    end
end

% savepdf(f, ['ICactivationStrength_AchillesLinear_' timeScale], '-dpdf')

end


function periodsFormat2 = changeFormat(periodsFormat1, behavior)

    periodsFormat2 = struct('PRE', [], 'RUN', [], 'POST', []);
    conditions = {'PRE', 'RUN', 'POST'};
    
    for ii = 1:3
       periodsFormat2.(conditions{ii}) = periodsFormat1(~(periodsFormat1(:,1) < behavior.time(ii,1) | periodsFormat1(:,2) > behavior.time(ii,2)), :); 
       periodsFormat2.(conditions{ii}) = periodsFormat2.(conditions{ii}) - behavior.time(ii,1);
    end
    
end


