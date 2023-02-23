function [paths, pathLens, pathPosRecon, jumpDisRL, maxJumpDisRL, jumpDisLR, maxJumpDisLR, posCoverRL, posCoverLR] = pathAnalysis(transmat, transProbThresh, overlapThresh, positionsRL, positionsLR, positionBins)

% finding all the paths(comprised of the states) based on the transition matrix (thresholded)
% originating from each node and removing all the redundant paths

[paths, pathLens] = pathFinder(transmat, transProbThresh, overlapThresh);

numofPaths = length(paths);

% for ip = 1 : length(paths)
%     
%     path = paths{ip};
%     
%     
%     % right 2 left
%     pathPosRecon{ip, 1} = positionsRL(:, path);
%     [~, maxPos] = max(pathPosRecon{ip, 1});
%     
%     maxPos = positionBins(maxPos);
%     currJumpsRL = abs(diff(maxPos)');
%     
%     jumpDisRL = [jumpDisRL; currJumpsRL];
%     maxJumpDisRL = [maxJumpDisRL max(currJumpsRL)];
%     
%     posCoverRL = [posCoverRL; max(maxPos) - min(maxPos)];
%     
%     
%     
%     % left 2 right
%     pathPosRecon{ip, 2} = positionsLR(:, path);
%     [~, maxPos] = max(pathPosRecon{ip, 2});
%     
%     maxPos = positionBins(maxPos);
%     currJumpsLR = abs(diff(maxPos)');
%     
%     jumpDisLR = [jumpDisLR; currJumpsLR];
%     maxJumpDisLR = [maxJumpDisLR max(currJumpsLR)];
%     
%     posCoverLR = [posCoverLR; max(maxPos) - min(maxPos)];
%     
% end
% 

% 
% 
% subplotdims = ceil(sqrt(numofPaths));
% figure;
% 
% for ip = 1 : numofPaths
%     subplot(subplotdims, subplotdims, ip)
%     imagesc(pathPosRecon{ip, 1});
% end
% 
% colormap(flipud(colormap('bone')))


end
