function plotICactivation4PBEs(eventsBinnedFiring, activationStrPBEs, posteriorProbMatrix, BDscore, ICscore, tempate_LR, template_RL, PBEidx2plot, fileBase, filename, textOnTop)

 binSize = 20; % ms
%%% set the configuration of the plots

plotheight = 50; %% in cm
plotwidth = 40;

nor = 4; %% number of events per row
noc = 10; %% number of events per column

leftedge = 1.2;    rightedge = 0.4;    topedge = 5;    bottomedge = 1.5;

gapr = 2;     gapc = 1;
fontsize = 5;



%%% set the positions of insets (rasters, bayesian posterior prob. matrix, and HMM likelihood profile) for each event

%%% each subplot (corresponding to an event) is comprised of 4 subplots

plotwidth2 = (plotwidth - (rightedge+leftedge) - (noc-1) * gapc)/noc;
plotheight2 = (plotheight - (topedge+bottomedge) - (nor-1) * gapr)/nor;

nor2 = 4;
noc2 = 1;

leftedge2 = 0;  rightedge2 = 0;     bottomedge2 = 0;    topedge2 = 0;

gapr2 = 0.2;    gapc2 = 1;



sub_pos = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); %% this gives us the positions for the subplots, each corresponded to an event


noPages = ceil(length(PBEidx2plot) / (nor * noc));
flag = 0; %% it is set to one when all the events are plotted


for plt = 1 : noPages
    
    if flag == 1 %%% do we need this part?
        break
    end
    
    
    %%% save and close the previous page
    
    if plt > 1
        
       filename=[fileBase '/ICsActivationStrengthDuringPBEs' num2str(plt-1)];
       saveas(gcf, [filename '.fig'])
       saveas(gcf, filename,  'epsc')
%        print(gcf, filename, '-dpng', '-r0')
       print(gcf, filename, '-dsvg', '-r0')

    end
    
%     close all
    
    
    %setting the Matlab figure for the current page
    
    f=figure;
    clf(f);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [plotwidth plotheight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
 
   
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

            
            currEvent  = eventsBinnedFiring{event, 1};
            currEvent  = currEvent(:, 1:floor(size(currEvent,2) / binSize) * binSize); %% truncate the events length to 
        
            spikesLR = currEvent(tempate_LR, :); %%% for raster of place cells in each direction
            spikesRL = currEvent(template_RL, :);
            


            sub_pos2 = subplot_pos_inset2(plotwidth2, plotheight2, leftedge2, rightedge2, bottomedge2, topedge2, noc2, nor2, gapr2, gapc2); %% this is a modified version of subplot_pos


            for ii = 1 : size(sub_pos2, 1)
                for jj = 1 : noc2

                    sub_pos2{ii,jj}(1) = (sub_pos2{ii,jj}(1) * plotwidth2 + sub_pos{i,j}(1) * plotwidth) / plotwidth;
                    sub_pos2{ii,jj}(2) = (sub_pos2{ii,jj}(2) * plotheight2 + sub_pos{i,j}(2) * plotheight) / plotheight;
                    sub_pos2{ii,jj}(3) = sub_pos2{ii,jj}(3) * plotwidth2 / plotwidth;
                    sub_pos2{ii,jj}(4) = sub_pos2{ii,jj}(4) * plotheight2 / plotheight;
                end
            end
            
           
            

            %%%% Spike raster plots
            
            
            %%% left to right template
            axes('position',sub_pos2{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

            for unit = 1 : size(spikesLR,1)
                
                spikes = find(spikesLR(unit,:));
                
                if ~isempty(spikes)
                    for spk = 1 : length(spikes)
                        plot([spikes(spk), spikes(spk)], [unit-1.25, unit+1.25], 'color', 'k', 'linewidth', 1) % colorset(ceil(tempate_LR(unit)*128/noUnits), :)
                        hold on
                    end
                end
                
            end
            
            xlim([0 size(spikesLR, 2)])
            ylim([0 length(tempate_LR)+3])
            
            if mod(eventInd, noc) == 1
               ylabel('Unit', 'fontsize', 6);
               set(gca, 'xtick',[], 'ytick',[1 size(spikesLR,1)], 'yticklabel', {'1',sprintf('%d', size(spikesLR,1))})
               title({['pbe ' num2str(event)]; 'leftright'},'fontsize', 6, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            
            else

                set(gca, 'xtick',[],'ytick',[])
                title({['pbe ' num2str(event)]; ''},'fontsize', 6, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            
            end
            
            %%% right to left template
            
            axes('position',sub_pos2{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

            for unit = 1 : size(spikesRL,1)
                
                spikes = find(spikesRL(unit,:));
                
                if ~isempty(spikes)
                    for spk = 1 : length(spikes)
                        plot([spikes(spk), spikes(spk)], [unit-1.25, unit+1.25], 'color', 'k', 'linewidth', 1) % colorset(ceil(tempate_LR(unit)*128/noUnits), :)
                        hold on
                    end
                end
                
            end
            
            xlim([0 size(spikesRL, 2)])
            ylim([0 length(template_RL)+3])
            
            set(gca, 'xtick',[],'ytick',[])
            
            hold on
            if mod(eventInd, noc) ~= 0  
                 
                 annotation('line',[sub_pos2{1,1}(1)+sub_pos2{1,1}(3)+0.005 sub_pos2{1,1}(1)+ sub_pos2{1,1}(3)+0.005], [sub_pos2{1,1}(2)+sub_pos2{4,1}(4) sub_pos2{4,1}(2)], 'linewidth', 2, 'color', [0.7 0.7 0.7])
            end   
            


            %%% Bayesian Decoding for each direction
            
            % LR direction
            
            ax1 = axes('position',sub_pos2{3,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
            
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
            
            hold on
            
            title(sprintf('\\color{magenta}UI-WC=%.1f\n\\color{blue}ICscore=%.0f', BDscore(event), ICscore(event)), 'fontsize', 5, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
     
            if mod(eventInd, noc) == 1
               ylabel('Position bin', 'fontsize', 6)
               set(ax1, 'xtick',[], 'ytick',[1 nPosBins], 'yticklabel', {'1',sprintf('%d', nPosBins)})
            else
               set(ax1, 'xtick', [], 'ytick', [])
            end
            
            
%             set(gca, 'visible', 'off')
            hold on
            

            %%% plot IC reactivation strength
            
            ax2 = axes('position',sub_pos2{4,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
            
%             nICs = size(activationStrPBEs{event}, 1);
%             plot(0.5:1:(size(activationStrPBEs{event}, 2)-0.5), activationStrPBEs{event}(setdiff(1:nICs, selectICs), :)', 'color', [0.7 0.7 0.7 0.3])
%             hold on
            plot(0.5:1:(size(activationStrPBEs{event}, 2)-0.5), activationStrPBEs{event}')
            
            
            xlim([0 size(activationStrPBEs{event}, 2)])
        end
    end
end


dim = [0.005 0.9 0.5 .1];

annotation('textbox',dim,'String',textOnTop, 'FitBoxToText','on'); %% , 'Interpreter', 'none'
filename = fullfile(fileBase, filename);
% filename=[fileBase '/pHMMnBD_tswap' num2str(plt)];

% filename = fileBase;
% saveas(gcf, [filename '.fig'])

saveas(gcf, filename,  'png')
saveas(gcf, filename,  'svg')
% saveas(gcf, filename,  'epsc')
print(gcf, filename, '-dpdf', '-r0')


end