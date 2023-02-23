function plotExamplePBEs2(eventBinnedfiring, sequenceScores,  statesProbDists, lspfPeakPositions, binsCongSignificance, PBEidx2plot,  binDur, fileinfo, behavior, FileBase)


binSize = binDur * 1000;
noUnits = size(eventBinnedfiring{1,1}, 1);


%%% set the configuration of the plots

plotheight = 40; %% in cm
plotwidth = 50;

nor = 3; %% number of events per row
noc = 10; %% number of events per column

leftedge = 1.2;    rightedge = 0.4;    topedge = 1;    bottomedge = 1.5;

gapr = 0; %0.5;
gapc = 0.5;
fontsize = 5;



%%% set the positions of insets (rasters, bayesian posterior prob. matrix, and HMM likelihood profile) for each event

%%% each subplot (corresponding to an event) is comprised of 4 subplots

plotwidth2 = (plotwidth - (rightedge+leftedge) - (noc-1) * gapc)/noc;
plotheight2 = (plotheight - (topedge+bottomedge) - (nor-1) * gapr)/nor;

nor2 = 13;
noc2 = 1;

leftedge2 = 0;  rightedge2 = 0;     bottomedge2 = 0;    topedge2 = 0;

gapr2 = 0.2;    gapc2 = 0.2;

colorset = colormap('jet'); %% colorset for the rasters
colorset = colorset(randperm(length(colorset)), :);


sub_pos = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); %% this gives us the positions for the subplots, each corresponded to an event


noPages = ceil(length(PBEidx2plot) / (nor * noc));
flag = 0; %% it is set to one when all the events are plotted

for plt = 1 : noPages
    
    if flag == 1 %%% do we need this part?
        break
    end
    
    
    %%% save and close the previous page
    
    if plt > 1
        
       filename=[FileBase '/examplePBEs_set' num2str(plt-1)];
       saveas(gcf, [filename '.fig'])
       print(gcf, filename, '-dpdf', '-r0')

    end
    
    close all
    
    
    %setting the Matlab figure for the current page
    
    f=figure;
    clf(f);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [plotwidth plotheight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
    
    
    
   
    pageFirstEventInd = (plt - 1) * 30 + 1; %%% index of the first event plotted in the current page
    
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
            

            event = PBEidx2plot(eventInd);

            
            currEvent = eventBinnedfiring{event, 1};
            currEvent = currEvent(:, 1:floor(size(currEvent,2) / binSize) * binSize); %% truncate the events length to 
        

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
            
            
            %%% leftward (RL) place cells
            axes('position',sub_pos2{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

            for unit = 1 : size(currEvent,1)
                
                spikes = find(currEvent(unit,:));
                
                if ~isempty(spikes)
                    for spk = 1 : length(spikes)
                        plot([spikes(spk), spikes(spk)], [unit-1, unit+1], 'color', colorset(ceil(unit*64/noUnits), :), 'linewidth', 1.5)
                        hold on
                    end
                end
                
            end
            
          
            xlim([0 size(currEvent, 2)])
            ylim([0 size(currEvent,1)+3])
            title(sprintf('PBE #%d, score %.2f', eventBinnedfiring{event, 3}, sequenceScores(event)), 'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            set(gca, 'xtick',[],'ytick',[])

            hold on

            

            %%% Bayesian Decoding for each direction
            ax1 = axes('position',sub_pos2{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
            
            
            %%% plot the posterior probability matrix
            imagesc(statesProbDists{event}, [0 1])
            set(gca,'YDir','normal')
            
            cm = flipud(bone);
            colormap(ax1, cm)
            set(gca,'ytick',[])
            
            hold on
            
            
            %%%
%             
%             axes('position',sub_pos2{3,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
%             
%             runIdx = find(fileinfo.xyt(:,3) > behavior.time(2,1) & fileinfo.xyt(:,3)< behavior.time(2,2));
%             runXpose = fileinfo.xyt(runIdx, 1);
%             runYpose = fileinfo.xyt(runIdx, 2);
%             
%             plot(runXpose-min(runXpose), runYpose-min(runYpose), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])
%             
%             
%             hold on
%             
%             PBEtrajectory = zeros(size(currEvent, 2)/20, 2);
%             for bin = 1:size(currEvent, 2)/20
%                 [~, maxState] = max(statesProbDists{event}(:, bin));
%                 PBEtrajectory(bin, :) = lspfPeakPositions(maxState, :);
%             end
%             
%             
%             binsDistThresh = 60;
%             maxnoSteps = plotTrajectory(PBEtrajectory, binsDistThresh);
%             
%             set(gca, 'xtick',[],'ytick',[])
%             title(sprintf('no steps = %d', maxnoSteps), 'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
%             
            
%             %%%% results of HMM (likelihood profiles)
%             ax2 = axes('position',sub_pos2{3,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
%             
% 
%             imagesc(binsCongSignificance{event}, [0 1])
%             colormap(ax2, 'bone')
%             set(gca, 'xtick',[],'ytick',[])


        end
    end
end

filename=[FileBase '/examplePBEs_set' num2str(plt)];%%% save the last page
saveas(gcf, [filename '.fig'])
print(gcf, filename, '-dpdf', '-r0')

end
