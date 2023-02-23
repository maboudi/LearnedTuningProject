function tempPlotBD3_1D(eventsBinnedFiring, posteriorProbMatrix, BDseqscore, weightedCorr, begPosition, endPosition, template, PBEidx2plot, binDur, fileBase, filename, textOnTop)



binSize = binDur * 1000;
% noUnits = size(eventsBinnedFiring{1,1}, 1);


%%% configure the plots

plotheight = 40; %% in cm
plotwidth  = 50;

nor = 3; %% number of events per row
noc = 20; %% number of events per column

leftedge = 1.2;    rightedge = 0.4;    topedge = 5;    bottomedge = 1.5;

gapr = 0.5;     gapc = 0.5;

fontsize = 5;


%%% set the positions of insets (rasters, bayesian posterior prob. matrix, and HMM likelihood profile) for each event

%%% each subplot (corresponding to an event) is comprised of 4 subplots

plotwidth2  = (plotwidth - (rightedge+leftedge) - (noc-1) * gapc)/noc;
plotheight2 = (plotheight - (topedge+bottomedge) - (nor-1) * gapr)/nor;

nor2 = 17;
% nor2 = 9;
noc2 = 1;

leftedge2 = 0;  rightedge2 = 0;     bottomedge2 = 0;    topedge2 = 0;

gapr2 = 0.2;    gapc2 = 0.2;


colorset = colormap('jet'); %% colorset for the rasters
colorset = colorset(randperm(length(colorset)), :);

colorset = [colorset; colorset];


sub_pos  = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); %% this gives us the positions for the subplots, each corresponded to an event


noPages = ceil(length(PBEidx2plot) / (nor * noc));
flag    = 0; %% it is set to one when all the events are plotted




for plt = 1 : noPages
    
    if flag == 1 %%% do we need this part?
        break
    end
    
    
    %%% save and close the previous page
    
    if plt > 1
        
       filename=[fileBase '/examplePBEs_set' num2str(plt-1)];
       saveas(gcf, [filename '.fig'])
       saveas(gcf, filename,  'epsc')
%        print(gcf, filename, '-dpng', '-r0')
       print(gcf, filename, '-dsvg', '-r0')

    end
    
    close all
    
    
    %setting the Matlab figure for the current page
    
    f=figure;
    clf(f);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [plotwidth plotheight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
    
%     
% 
%     for ii = 1: nor
%         for jj = 1: noc-1%
%             line([leftedge + jj*plotwidth/noc leftedge + jj*plotwidth/noc], [topedge + (ii-1)*plotheight/nor+3 topedge + ii*plotheight/nor-3], 'linewidth', 2, 'color', [0.7 0.7 0.7])
% 
%         end
%     end
    
   
    pageFirstEventInd = (plt - 1) * nor*noc + 1; %%% index of the first event plotted in the current page
    
    for i = 1:nor
        
        if flag == 1 %%% leave the loop when all the events were plotted
            break
        end
        
        for j = 1:noc

            eventInd = (i - 1) * noc + j + pageFirstEventInd - 1;
            
            
            if eventInd > length(PBEidx2plot) %% leave if all the events were plotted
                flag = 1;
                break
            end
            

            event     = PBEidx2plot(eventInd); 
            eventID   = eventsBinnedFiring{event, 3};
            
            currEvent  = eventsBinnedFiring{event, 1};
            currEvent  = currEvent(:, 1:floor(size(currEvent,2) / binSize) * binSize); %% truncate the events length to 
        
%             spikesL2R = currEvent(L2R_template, :); %%% for raster of place cells in each direction
            spikes = currEvent(template, :);
            


            sub_pos2 = subplot_pos2(plotwidth2, plotheight2, leftedge2, rightedge2, bottomedge2, topedge2, noc2, nor2, gapr2, gapc2); %% this is a modified version of subplot_pos


            for ii = 1 : size(sub_pos2, 1)
                for jj = 1 : noc2

                    sub_pos2{ii,jj}(1) = (sub_pos2{ii,jj}(1) * plotwidth2 + sub_pos{i,j}(1) * plotwidth) / plotwidth;
                    sub_pos2{ii,jj}(2) = (sub_pos2{ii,jj}(2) * plotheight2 + sub_pos{i,j}(2) * plotheight) / plotheight;
                    sub_pos2{ii,jj}(3) = sub_pos2{ii,jj}(3) * plotwidth2 / plotwidth;
                    sub_pos2{ii,jj}(4) = sub_pos2{ii,jj}(4) * plotheight2 / plotheight;
                end
            end
            
           
            

            %%%% Spike raster plots
            
            
            axes('position',sub_pos2{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

            for unit = 1 : size(spikes,1)
                
                currspikes = find(spikes(unit,:));
                
                if ~isempty(currspikes)
                    for spk = 1 : length(currspikes)
                        plot([currspikes(spk), currspikes(spk)], [unit-1.25, unit+1.25], 'color', 'k', 'linewidth', 1) % colorset(ceil(L2R_template(unit)*128/noUnits), :)
                        hold on
                    end
                end
                
            end
            
            xlim([0 size(spikes, 2)])
            ylim([0 length(template)+3])
            
            if mod(eventInd, noc) == 1
               ylabel('Unit', 'fontsize', 6);
               set(gca, 'xtick',[], 'ytick',[1 size(spikes,1)], 'yticklabel', {'1',sprintf('%d', size(spikes,1))})
               title({['pbe ' num2str(eventID)], 'Template 1'},'fontsize', 6, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            
            else
                
                set(gca, 'xtick',[],'ytick',[])
                title({['pbe ' num2str(eventID)], ''},'fontsize', 6, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            
            end
            hold on
            if mod(eventInd, noc) ~= 0  
                 
                 annotation('line',[sub_pos2{1,1}(1)+sub_pos2{1,1}(3)+0.005 sub_pos2{1,1}(1)+ sub_pos2{1,1}(3)+0.005], [sub_pos2{1,1}(2)+sub_pos2{3,1}(4) sub_pos2{3,1}(2)], 'linewidth', 2, 'color', [0.7 0.7 0.7])
            end   

            
            
            %%% Bayesian Decoding

            
            ax1 = axes('position',sub_pos2{3,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
%             ax1 = axes('position',sub_pos2{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
            
            currPPM   = posteriorProbMatrix{event,1};
            nPosBins = size(currPPM, 1);
            
            maxProb   = max(currPPM(:));
            minProb   = min(currPPM(:));
            
            Clim_low  = minProb;
            Clim_high = minProb + 0.75 *(maxProb - minProb);
            
            
            %%% plot the posterior probability matrix
            imagesc(currPPM, [Clim_low Clim_high])
            set(gca,'YDir','normal')
            colormap(ax1, 'hot')
            title({sprintf('TS-RT=%.1f\nTS-WC=%.1f(%.2f)', BDseqscore.wPBEtimeswap.replayScore.prctilescore(event, 1) , BDseqscore.wPBEtimeswap.weightedCorr.prctilescore(event, 1), weightedCorr(event, 1)); ... 
                sprintf('UI-RT=%.1f\nUI-WC=%.1f', BDseqscore.unitIDshuffle.replayScore.prctilescore(event, 1) , BDseqscore.unitIDshuffle.weightedCorr.prctilescore(event, 1))}, 'fontsize', 5, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            
            
            hold on
            
        
            %%%% plot the best fitting line
%             
%             plot([begPosition(event,1,1) endPosition(event,1,1)], [begPosition(event,2,1) endPosition(event,2,1)],'color','w', 'linewidth', 0.5)
%             
%             plot([begPosition(event,1,1) endPosition(event,1,1)], [begPosition(event,2,1)-12 endPosition(event,2,1)-12],'color','w', 'linestyle', ':', 'linewidth', 0.5)
%             plot([begPosition(event,1,1) endPosition(event,1,1)], [begPosition(event,2,1)+12 endPosition(event,2,1)+12],'color','w', 'linestyle', ':', 'linewidth', 0.5)
%             
            if mod(eventInd, noc) == 1
               ylabel('Position bin', 'fontsize', 6)
               set(ax1, 'xtick',[], 'ytick',[1 nPosBins], 'yticklabel', {'1',sprintf('%d', nPosBins)})
            else
               set(ax1, 'xtick', [], 'ytick', [])
            end
            
%             set(gca, 'visible', 'off')
            hold on
  

        end
    end
end


dim = [0.005 0.85 0.5 .1];

annotation('textbox',dim,'String',textOnTop, 'FitBoxToText','on', 'Interpreter', 'none');

filename = fullfile(fileBase, filename);
% filename = [fileBase '/pHMMnBD_tswap' num2str(plt)];

% filename = fileBase;
% saveas(gcf, [filename '.fig'])
% saveas(gcf, filename,  'epsc')

print(gcf, filename, '-dpdf', '-r0')

% close all
end


