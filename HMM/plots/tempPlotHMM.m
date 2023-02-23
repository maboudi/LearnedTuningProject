function tempPlotHMM(eventsBinnedFiring, statesProbDists, seqScore, unitsSortIdx, PBEidx2plot, binDur, fileBase)



binSize = binDur * 1000;
% noUnits = size(eventsBinnedFiring{1,1}, 1);


%%% configure the plots

plotheight = 40; %% in cm
plotwidth  = 50;

nor = 5; %% number of events per row
noc = 20; %% number of events per column

leftedge = 1.2;    rightedge = 0.4;    topedge = 1;    bottomedge = 1.5;

gapr = 0.5;     gapc = 0.5;

fontsize = 5;


%%% set the positions of insets (rasters, bayesian posterior prob. matrix, and HMM likelihood profile) for each event

%%% each subplot (corresponding to an event) is comprised of 4 subplots

plotwidth2  = (plotwidth - (rightedge+leftedge) - (noc-1) * gapc)/noc;
plotheight2 = (plotheight - (topedge+bottomedge) - (nor-1) * gapr)/nor;

% nor2 = 17;
nor2 = 9;
noc2 = 1;

leftedge2 = 0;  rightedge2 = 0;     bottomedge2 = 0;    topedge2 = 0;

gapr2 = 0.2;    gapc2 = 0.2;


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
            

            event = PBEidx2plot(eventInd);

            
            currEvent  = eventsBinnedFiring{event, 1};
            currEvent  = currEvent(:, 1:floor(size(currEvent,2) / binSize) * binSize); 
            
            
%             spikes = currEvent(unitsSortIdx(:, event), :); %%% for raster of place cells in each direction
%             spikesR2L = currEvent(R2L_template, :);
            spikes = currEvent(unitsSortIdx, :);
            
            sub_pos2 = subplot_pos3(plotwidth2, plotheight2, leftedge2, rightedge2, bottomedge2, topedge2, noc2, nor2, gapr2, gapc2); %% this is a modified version of subplot_pos


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

            for unit = 1:size(spikes,1)
                
                currspikes = find(spikes(unit,:));
                if ~isempty(currspikes)
                    for spk = 1 : length(currspikes)
                        plot([currspikes(spk) currspikes(spk)], [unit-1.25 unit+1.25], 'color', 'k','linewidth', 1) %%, colorset(ceil(R2L_template(unit)*128/noUnits), :)
                        hold on
                    end
                end
                
            end
            
            title(sprintf('sequence score %.1f', seqScore(event)),'fontsize', 5, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            
            xlim([0 size(spikes, 2)])
%             ylim([0 length(template)+3])
            set(gca, 'xtick',[],'ytick',[])
            hold on
            
            
            
            %%% states probability distributions for each PBE
                   
            
            ax2 = axes('position',sub_pos2{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
            
            currProbDist   = statesProbDists{event};
            
            maxProb   = max(currProbDist(:));
            minProb   = min(currProbDist(:));
            
            Clim_low  = minProb;
            Clim_high = minProb + 0.5 *(maxProb - minProb);
            
            %%% plot the posterior probability matrix
            imagesc(currProbDist, [Clim_low Clim_high])
            set(gca,'YDir','normal')
            colormap(ax2, 'jet')
            
            hold on

            set(gca, 'xtick',[],'ytick',[])

            
            
        end
    end
end
            

filename=[fileBase '/examplePBEs_set' num2str(plt)];
saveas(gcf, [filename '.fig'])
saveas(gcf, filename,  'epsc')
print(gcf, filename, '-dsvg', '-r0')

% close all
end


