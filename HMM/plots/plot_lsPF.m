function lsPFs = plot_lsPF(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod)


% numofStates = 40;
lsPFs = HMM_StatesPlaceField_2D(runBinnedfiring{2}, transmat, lambda, posbinIdx, numofStates, xposcenters, yposcenters);



% plot the lsPFs

figure;

numStatesperRow = 3;

runidx = find(fileinfo.xyt(:, 3) > behavior.time(2,1) & fileinfo.xyt(:, 3) < behavior.time(2,2));
xpos = fileinfo.xyt(runidx,1);
ypos = fileinfo.xyt(runidx,2);

for ii = 1: numofStates
    
    subplot(ceil(numofStates/numStatesperRow), numStatesperRow, ii)
    
    plot(fileinfo.xyt(runidx, 1), fileinfo.xyt(runidx, 2), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])
    
    hold on
    h1 = imagesc(xposcenters, yposcenters, lsPFs(:,:, ii));
    set(h1, 'AlphaData', 0.7)

    
    colormap('hot')
    set(gca,'YDir','normal')
    
    cornerstate = (ceil(numofStates/numStatesperRow)-1) * numStatesperRow + 1;
    
    if ii == cornerstate
        
        set(gca, 'XTick', floor([min(xpos) max(xpos)]), 'YTick', floor([min(ypos) max(ypos)]), 'box', 'off')
        xlabel('X (cm)', 'fontsize', 14)
        ylabel('Y (cm)', 'fontsize', 14)
        
    else
        
        set(gca, 'XTick', [], 'YTick', [], 'box', 'off')
    end

    xlim(floor([min(xpos) max(xpos)]))
    ylim(floor([min(ypos) max(ypos)]))

end

% suptitle([fileinfo.name '-' pbePeriod])
set(gcf, 'position', [1 1 800 1500])

end
