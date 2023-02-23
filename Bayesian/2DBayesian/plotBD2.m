function plotBD2(eventBinnedfiring, peakPosBin, distFromPrevBin, PBETrajsegments, PBEidx2plot, fileinfo, behavior, binDur, posBinSize, FileBase)


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

nor2 = 9;
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
%        saveas(gcf, [filename '.fig'])
       print(gcf, filename, '-dpdf', '-r0')
       print(gcf, filename, '-dpng', '-r0')

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
        

            sub_pos2 = subplot_pos3(plotwidth2, plotheight2, leftedge2, rightedge2, bottomedge2, topedge2, noc2, nor2, gapr2, gapc2); %% this is a modified version of subplot_pos
            
            
            
            for ii = 1 : size(sub_pos2, 1)
                for jj = 1 : noc2

                    sub_pos2{ii,jj}(1) = (sub_pos2{ii,jj}(1) * plotwidth2 + sub_pos{i,j}(1) * plotwidth) / plotwidth;
                    sub_pos2{ii,jj}(2) = (sub_pos2{ii,jj}(2) * plotheight2 + sub_pos{i,j}(2) * plotheight) / plotheight;
                    sub_pos2{ii,jj}(3) = sub_pos2{ii,jj}(3) * plotwidth2 / plotwidth;
                    sub_pos2{ii,jj}(4) = sub_pos2{ii,jj}(4) * plotheight2 / plotheight;
                end
            end
            
            
            
            % trajectories to plot
            
            currPBEtrajs = PBETrajsegments{event};

            currPositions = peakPosBin{event}*posBinSize;
            currPositions = currPositions(:,2:-1:1);

            currDistances = distFromPrevBin{event};

            noTimeBins = size(currPositions, 1);

            trajs2Plot = find(currPBEtrajs(:,3) > 2); % both criteria on minimum number of steps and distance covered were met
%  & currPBEtrajs(:,5) > 60

            

%             hold on

            
            
            % plot the virtual trajectories
        
            axes('position',sub_pos2{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
            
            runIdx = find(fileinfo.xyt(:,3) > behavior.time(2,1) & fileinfo.xyt(:,3)< behavior.time(2,2));
            runXpose = fileinfo.xyt(runIdx, 1);
            runYpose = fileinfo.xyt(runIdx, 2);
            
            plot(runXpose-min(runXpose), runYpose-min(runYpose), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])
            
            
            hold on
            
            
            % plot overlaid trajectories
            
           colorset = colormap('jet');
        
            for tt = 1 : length(trajs2Plot)
                currTraj = currPBEtrajs(trajs2Plot(tt), 1:2);
                noSteps = currPBEtrajs(trajs2Plot(tt), 3);

                bins2plot = currDistances(find(currDistances(:, 2) == currTraj(1)):find(currDistances(:, 2) == currTraj(2)), 2);

                for ii = 1: noSteps+1

                    if ii > 1

                        currBin = bins2plot(ii);
                        preBin = bins2plot(ii-1);
                       line([currPositions(currBin, 1) currPositions(preBin, 1)], [currPositions(currBin, 2) currPositions(preBin, 2)], 'color', colorset(ceil(currBin/noTimeBins*size(colorset, 1)), :), 'linewidth', 2)

                        plot(currPositions(currBin, 1), currPositions(currBin, 2), '.', 'markersize', 15, 'color', colorset(ceil(currBin/noTimeBins*size(colorset, 1)), :))
                        plot(currPositions(preBin, 1), currPositions(preBin, 2), '.', 'markersize', 15, 'color', colorset(ceil(preBin/noTimeBins*size(colorset, 1)), :))

                        
                    end
                end
                
                text(currPositions(bins2plot(1), 1)-10, currPositions(bins2plot(1), 2)+5, sprintf('%d', bins2plot(1)), 'fontsize', 7, 'color', 'b');

                text(currPositions(bins2plot(end), 1)-10, currPositions(bins2plot(end), 2)+5, sprintf('%d', bins2plot(end)), 'fontsize', 7, 'color', 'k');


            end

            set(gca, 'xtick',[],'ytick',[])
            

            xlim([0 max(runXpose-min(runXpose))])
            ylim([0 max(runYpose-min(runYpose))])
            
            
            hold on
            
            %%%% Spike raster plots
            
            colorset = colorset(randperm(length(colorset)), :);

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
            title(sprintf('PBE #%d, no steps=%d', eventBinnedfiring{event, 3}, max(currPBEtrajs(:, 3))), 'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
            
            
            for tt = 1 : length(trajs2Plot)
                currTraj = currPBEtrajs(trajs2Plot(tt), 1:2);
                noSteps = currPBEtrajs(trajs2Plot(tt), 3);
                
                bins2plot = currDistances(find(currDistances(:, 2) == currTraj(1)):find(currDistances(:, 2) == currTraj(2)), 2);
                
                hold on
                line([(bins2plot(1)-1)*binSize+1 (bins2plot(1)-1)*binSize+1], [0 size(currEvent,1)+3], 'color', 'b')
                
                hold on
                line([(bins2plot(end)-1)*binSize+1 (bins2plot(end)-1)*binSize+1], [0 size(currEvent,1)+3], 'color', 'k')
                
            end
                

            
            set(gca, 'xtick',[],'ytick',[])
            
        end
    end
end

filename=[FileBase '/examplePBEs_set' num2str(plt)];%%% save the last page
% saveas(gcf, [filename '.fig'])
print(gcf, filename, '-dpdf', '-r0')
print(gcf, filename, '-dpng', '-r0')
end
            
            
            
    

        