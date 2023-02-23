function maxnoSteps = plotTrajectory(positions, frame2framedisThresh)

colorset = flipud(jet);

X = positions;
noBins = size(X, 1);

% flag = 0;
% maxnoSteps = 0;

distFromPrevBin = zeros(noBins, 1);
for ii = 1:noBins
    
    if ii > 1
        distFromPrevBin(ii) = norm(X(ii, :) - X(ii-1, :));
    end
    
end

distFromPrevBin = [distFromPrevBin (1:noBins)'];

distFromPrevBin(distFromPrevBin(:,1) < 1, :) = []; % removing all the stationary segments
distFromPrevBin = [0 1; distFromPrevBin];

isDistLessThresh = distFromPrevBin(2:end, 1) <= frame2framedisThresh; % in cm
isDistLessThresh = [0; isDistLessThresh]; 

% find the consequitive time bins with bin t0 bin distance less than
% the threshold

temp = abs(diff([0; isDistLessThresh])); 

startOrEndBins = find(temp > 0)-1; % decrement by 1 
    
if mod(length(startOrEndBins), 2) == 1
    startOrEndBins = [startOrEndBins; size(distFromPrevBin, 1)];
end

numSegs = length(startOrEndBins)/2;

if numSegs == 0
    maxnoSteps = 0;
    return
end

PBETrajsegments = zeros(numSegs, 4); % start and end of each (potentially) trajectory segments, without filtering based on the length yet
for seg = 1: numSegs

    segStart  = startOrEndBins((seg-1)*2+1);
    segEnd    = startOrEndBins(seg*2);
    segNSteps = segEnd - segStart;

    segDistanceCovered = norm(X(distFromPrevBin(segEnd, 2), :) - X(distFromPrevBin(segStart, 2), :));
    
    
    PBETrajsegments(seg, :) = [distFromPrevBin(segStart, 2) distFromPrevBin(segEnd, 2) ...
                                    segNSteps  ... % trajectory's number of steps
                                    segDistanceCovered]; % trajectory's covered distance


end

PBETrajsegments(PBETrajsegments(:, 4) < 1, :) = [];

if isempty(PBETrajsegments)
    maxnoSteps = 0;
    return
end

longestTrajInd = find(PBETrajsegments(:, 3) == max(PBETrajsegments(:, 3)));

if numel(longestTrajInd) > 1
    
    distance = PBETrajsegments(longestTrajInd, 4);

    longestTrajInd = longestTrajInd(distance == max(distance));
    
end


trajectory2Plot = PBETrajsegments(longestTrajInd, 1:2);
maxnoSteps = PBETrajsegments(longestTrajInd, 3);

bins2plot = distFromPrevBin(find(distFromPrevBin(:, 2) == trajectory2Plot(1)):find(distFromPrevBin(:, 2) == trajectory2Plot(2)), 2);

for ii = 1: maxnoSteps+1
    
    if ii > 1
        
        currBin = bins2plot(ii);
        preBin = bins2plot(ii-1);
       line([X(currBin, 1) X(preBin, 1)], [X(currBin, 2) X(preBin, 2)], 'color', colorset(ceil(currBin/noBins*size(colorset, 1)), :), 'linewidth', 2)
        
        plot(X(currBin, 1), X(currBin, 2), '.', 'markersize', 15, 'color', colorset(ceil(currBin/noBins*size(colorset, 1)), :))
        plot(X(preBin, 1), X(preBin, 2), '.', 'markersize', 15, 'color', colorset(ceil(preBin/noBins*size(colorset, 1)), :))
        
    end
    
end

text(X(bins2plot(1), 1)-10, X(bins2plot(1), 2)+5, sprintf('%d', bins2plot(1)), 'fontsize', 7, 'color', 'b');

text(X(bins2plot(end), 1)-10, X(bins2plot(end), 2)+5, sprintf('%d', bins2plot(end)), 'fontsize', 7, 'color', 'k');



end


