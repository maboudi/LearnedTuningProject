function pathsPlot2(paths, pathPosRecon, positionBins)


numofPaths = length(paths);

subplotdims = ceil(sqrt(numofPaths));

directions = {'R2L Direction','L2R Direction'};

myColors = colormap('jet');
myColors(1,:) = [1 1 1];


    
for dir = 1:2 % RL and LR
    
    figure;
    
    for ip =  1 : numofPaths
        
        subplot(subplotdims, subplotdims, ip)
        
        curr_path = paths{ip};
        curr_spaceMap = pathPosRecon{ip, dir};  
        
        
        hold on
        
        for istate = 1 :  length(curr_path)
            
            curr_state = curr_path(istate);
            
            curr_map = curr_spaceMap(:, istate);
            curr_map = curr_map./max(curr_map);
            
            [maxP, idx] = max(curr_map);
            
            h = fill([curr_map+3*istate; flipud(3*istate - curr_map)], [positionBins' ;flipud(positionBins')], myColors(curr_state, :)); 
            set(h, 'EdgeColor', 'none')
            
            plot([maxP+3*istate 3*istate-maxP], [positionBins(idx) positionBins(idx)], 'y','linewidth', 1)
            
        end
        
        hold off
        
        set(gca, 'XTick', [], 'YTick', [])
        xlim([0 3*15])
        xlabel(num2str(ip), 'fontsize', 9)
        
    end
        
    annotation('textbox', [0 0.9 1 0.1], ...
        'String', directions{dir}, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'fontsize', 28)

end


end


