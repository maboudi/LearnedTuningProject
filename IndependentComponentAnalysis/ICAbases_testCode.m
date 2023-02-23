
PRE_data = cell2mat(eventsBinnedfiringPRE(:, 2)');
RUN_data = cell2mat(eventsBinnedfiringRUN(:, 2)');
POST_data = cell2mat(eventsBinnedfiringPOST(:, 2)');

% 
% for ii = 1: size(PRE_data, 1)
%     PRE_data(ii,:) = PRE_data(ii, randperm(size(PRE_data, 2)));
% end
% 
% for ii = 1: size(RUN_data, 1)
%     RUN_data(ii,:) = RUN_data(ii, randperm(size(RUN_data, 2)));
% end
% 
% for ii = 1: size(POST_data, 1)
%     POST_data(ii,:) = POST_data(ii, randperm(size(POST_data, 2)));
% end

PRE_ICs = pca_ica(PRE_data, 10);
RUN_ICs = pca_ica(RUN_data, 10);
POST_ICs = pca_ica(POST_data, 10);


% PRE and RUN

[RUN_PRE_PopVecSim, RUN_PRE_lambdaSortIdx] = PopVecSimilarity(RUN_ICs, PRE_ICs);


% RUN and POST

[RUN_POST_PopVecSim, RUN_POST_lambdaSortIdx] = PopVecSimilarity(RUN_ICs, POST_ICs);


% PRE and POST

[PRE_POST_PopVecSim, PRE_POST_lambdaSortIdx] = PopVecSimilarity(PRE_ICs, POST_ICs);



figure; 
set(gcf, 'position', [-1528 408 1185 268])

subplot(1,3,1)
imagesc(RUN_PRE_PopVecSim(RUN_PRE_lambdaSortIdx, :), [40 90])
title('RUN and PRE')

subplot(1,3,2)
imagesc(RUN_POST_PopVecSim(RUN_POST_lambdaSortIdx, :), [40 90])
title('RUN and POST')

subplot(1,3,3)
imagesc(PRE_POST_PopVecSim(PRE_POST_lambdaSortIdx, :), [40 90])
title('PRE and POST')
origSize = get(gca, 'position');

colorbar
set(gca, 'position', origSize)

mymap = colormap('copper');
mymap = flipud(mymap);
colormap(mymap)



function [PopVectSim, sortAstates] = PopVecSimilarity(x, y)

PopVectSim = acos(x' * y ./ (repmat(diag(sqrt(x'*x)), [1 size(y,2)]) .*  repmat(diag(sqrt(y'*y))', [size(x,2) 1])))*180/pi;

PopVectSim2 = 180 - PopVectSim;   
PopVectSim(PopVectSim > 90) = PopVectSim2(PopVectSim > 90);

[~, maxBstates] = min(PopVectSim, [], 2);
[~, sortAstates] = sort(maxBstates);

% PopVectSim = PopVectSim(sortAstates, :);

end


