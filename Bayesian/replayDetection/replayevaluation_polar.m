function [replayScore, weightedCorr, jumpDistance, begPnt, endPnt] = replayevaluation_polar(postPr)


% weighted correlation
% weightedCorr = [];


weightedCorr = calWeightedCorr(postPr);


% find the best fit line

% tic

[nPosBins, nTimeBins] = size(postPr);

t_mid = (nTimeBins - 1)/2;
p_mid = (nPosBins - 1)/2;

diam = sqrt((nTimeBins-1)^2 + (nPosBins-1)^2);

minAngle = pi/2 - atan(nPosBins/4);

thetas      = linspace(minAngle, 2*pi-minAngle, 100);
sin_thetas  = sin(thetas);
cos_thetas  = cos(thetas);

% thetas  = linspace(-pi/2, pi/2, 150);
lengths = linspace(0, diam/3, 50);



% replacing each element within a column with a summation over the neighboring elements
% (corresponding to 30 sm vicinity (-15 cm to 15 cm))

nVicinityBins = 7; 
postPr_summed = zeros(size(postPr));
for ii = 1: nTimeBins
    postPr_summed(:, ii) = conv(postPr(:, ii), ones(2*nVicinityBins+1, 1), 'same');
end


sumProb = zeros(numel(thetas), numel(lengths));


for ii = 1: numel(thetas)
    for jj = 1: numel(lengths)

        tempSum = 0;
        for t = 1:nTimeBins
            
            p = ceil((lengths(jj) - (t-t_mid)*cos_thetas(ii))/sin_thetas(ii) + p_mid);
            if p <= nPosBins && p >= 1
               tempSum = tempSum + postPr_summed(p, t);
            else
               tempSum = tempSum + median(postPr_summed(:, t));
            end
            
        end
        
        sumProb(ii, jj) = tempSum/nTimeBins; 

    
    end
end



%%% finding the slope and rho (y intercept) resulting in the highest
%%% goodness of fit

replayScore = max(sumProb(:));

[ind1, ind2] = find(sumProb == replayScore, 1, 'first');
bestTheta   = thetas(ind1);
bestLength  = lengths(ind2);


%%% start and end points of the best fitting line

tBeg = 1;
yBeg = floor((bestLength - (tBeg-t_mid)*cos(bestTheta))/sin(bestTheta) + p_mid);

begPnt = [tBeg yBeg];



tEnd = nTimeBins;
yEnd = floor((bestLength - (nTimeBins-t_mid)*cos(bestTheta))/sin(bestTheta) + p_mid);

endPnt = [tEnd yEnd];
% toc
% 
% slope = tan(bestTheta);
% rho   =  yBeg ;

% 
% 
% figure; imagesc(postPr); colormap('jet'); set(gca, 'yDir', 'normal')
% hold on
% line([begPnt(1) endPnt(1)], [begPnt(2) endPnt(2)], 'linewidth', 3, 'color', 'w')
% hold on
% text(t_mid, p_mid, '*', 'fontsize', 10, 'color', 'w')
% 
% figure; imagesc(sumProb); colormap('jet'); set(gca, 'yDir', 'normal')

%% jump distance
% tic
[~, peakPos] = max(postPr, [], 1);


% eliminating non-firing time bins
sumProb                = sum(postPr, 1);
peakPos(sumProb < 0.1) = [];


% calculating jump distances



if numel(peakPos) < 2
   
   maxJump     = 1;
   medianJump  = 1;
   normMaxJump = 1;
   
else

    jumpDistances = abs(diff(peakPos));
    
    maxJump       = max(jumpDistances)/nPosBins; % maximum jump distance as a ratio of the track length
    medianJump    = median(jumpDistances)/nPosBins; % median jump distance as a ratio of the track length


    % normalized maximum jump distance (data as a percentile of shuffle distribution:within-PBE time bin shuffle (ref: Dragoi2019 Science))

    nShuffles      = 500; 
    shuffleMaxJump = zeros(nShuffles, 1);

    for ii = 1: nShuffles

        shufflePeakPos      = peakPos(randperm(length(peakPos)));
        shuffleJumpDistance = abs(diff(shufflePeakPos));


        shuffleMaxJump(ii)  = max(shuffleJumpDistance)/nPosBins;

    end

    normMaxJump  = length(find(shuffleMaxJump < maxJump))/nShuffles;

end


jumpDistance = [maxJump normMaxJump medianJump];

% toc


end