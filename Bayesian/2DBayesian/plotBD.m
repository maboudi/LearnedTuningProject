function plotBD(peakPosBin, distFromPrevBin, PBETrajsegments, fileinfo, behavior, posBinSize)

% this function is used for plotting 2D Bayesian trajectories

%%% inputs:

% peakPosBin: peak position bins corresponding to each 20 ms time bin
% distFromPrevBin: difference in peak position bin from the previous time
% bin


figure;
nPBEs2Plot = 64;
colorset = flipud(jet);

nPBEs = length(PBETrajsegments);

t = 0; % indices of PBEs to plot

pbeind = randperm(nPBEs);

for ii = 1: nPBEs
    
    pbe = pbeind(ii);
    
    if t > nPBEs2Plot
        break
    end
    
    currPBEtrajs = PBETrajsegments{pbe};
    
    currPositions = peakPosBin{pbe}*posBinSize;
    currPositions = currPositions(:,2:-1:1);
    
    currDistances = distFromPrevBin{pbe};
    
    noTimeBins = size(currPositions, 1);
    
    trajs2Plot = find(currPBEtrajs(:,3) > 2 & currPBEtrajs(:,5) > 60); % both criteria on minimum number of steps and distance covered were met
    
    if ~isempty(trajs2Plot)
        
        t = t + 1; 
        
        if t > nPBEs2Plot
            break
        end
        
        subplot(8,8,t)

        runIdx = find(fileinfo.xyt(:,3) > behavior.time(2,1) & fileinfo.xyt(:,3)< behavior.time(2,2));
        runXpose = fileinfo.xyt(runIdx, 1);
        runYpose = fileinfo.xyt(runIdx, 2);

        plot(runXpose-min(runXpose)+posBinSize/2, runYpose-min(runYpose)+posBinSize/2, '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])

        hold on
        
        for tt = 1 : length(trajs2Plot)
            currTraj = currPBEtrajs(trajs2Plot(tt), 1:2);
            noSteps = currPBEtrajs(trajs2Plot(tt), 3);
            
            bins2plot = currDistances(find(currDistances(:, 2) == currTraj(1)):find(currDistances(:, 2) == currTraj(2)), 2);
            
            for jj = 1: noSteps+1
                
                if jj > 1
                   
                    currBin = bins2plot(jj);
                    preBin = bins2plot(jj-1);
                   line([currPositions(currBin, 1) currPositions(preBin, 1)], [currPositions(currBin, 2) currPositions(preBin, 2)], 'color', colorset(ceil(currBin/noTimeBins*size(colorset, 1)), :), 'linewidth', 2)

                    plot(currPositions(currBin, 1), currPositions(currBin, 2), '.', 'markersize', 15, 'color', colorset(ceil(currBin/noTimeBins*size(colorset, 1)), :))
                    plot(currPositions(preBin, 1), currPositions(preBin, 2), '.', 'markersize', 15, 'color', colorset(ceil(preBin/noTimeBins*size(colorset, 1)), :))
                    
                    
                end
            end  
        end
        
        set(gca, 'xtick',[],'ytick',[])
        title(sprintf('PBE # %d, no steps = %d', pbe, max(currPBEtrajs(:, 3))), 'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
        
        xlim([0 max(runXpose-min(runXpose))])
        ylim([0 max(runYpose-min(runYpose))])
    end

end

end
