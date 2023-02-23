function plotPanel(stateCorrespMat, rowSortInd, model, model2comp, ccodeRange, currTitle)


% 
% for t = 1:5
%     [i,j] = find(stateCorrespMat == max(stateCorrespMat(:)));
%     stateCorrespMat(i,j) = mean(stateCorrespMat(:));
% end

if ~isempty(ccodeRange)
    
    imagesc(stateCorrespMat(rowSortInd, :), ccodeRange); colormap('jet'); % 
else
    imagesc(stateCorrespMat(rowSortInd, :)); colormap('jet'); % 
end


set(gca, 'YDir', 'normal', 'ytick', 1:length(rowSortInd), 'yticklabel', rowSortInd, 'fontsize', 10)
set(gca, 'XTick', [], 'YTick', [])

ylabel([model ' states'] , 'fontsize', 10)
xlabel([model2comp ' states'], 'fontsize', 10)

title(currTitle, 'fontsize', 10, 'fontweight', 'normal')
% colorbar

end
