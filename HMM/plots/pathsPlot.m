function pathsPlot(paths, pathPosRecon, positionBins)


numofPaths = length(paths);

subplotdims = ceil(sqrt(numofPaths));

directions = {'R2L Direction','L2R Direction'};
dirColors = {'r','b'};

myColors = colormap('jet');
myColors(1,:) = [1 1 1];
    
    
for dir = 1:2 % RL and LR
    
    figure;
    
    for ip =  1 : numofPaths

        curr_path = paths{ip};

        curr_spaceMap = pathPosRecon{ip, dir};   
        curr_spaceMap(find(curr_spaceMap < 0.025)) = 0;

        curr_spaceMap = curr_spaceMap + 10*repmat(curr_path, size(curr_spaceMap, 1), 1); 
        curr_spaceMap(find(floor(curr_spaceMap) == curr_spaceMap)) = 0;
        
        
        subplot(subplotdims, subplotdims, ip)
        imagesc(1:length(curr_path), positionBins, curr_spaceMap); set(gca,'YDir','normal')
        set(gca, 'XTick', [], 'YTick', [])
        
        if ip == 1

           text(0, 0, 'Position', 'rotation', 90, 'fontsize', 24)
           text(1, 300, 'State', 'fontsize', 24)
        end

    end
    
    colormap(myColors)
    annotation('textbox', [0 0.9 1 0.1], ...
        'String', directions{dir}, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'fontsize', 28)

end


end


