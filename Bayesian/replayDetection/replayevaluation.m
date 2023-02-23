function [maxgoodnessofFit, bestSlope, bestrho, begPnt, endPnt] = replayevaluation(postPr, varargin)

% This function is intended to find the best (virtual) trajectory based on
% the posterior probability matrix. To do so, a greedy search is done
% to find the line with the maximum goodness of fit to the matrix. 

[noPosBins, noTimeBins]  = size(postPr);


% Define the search ranges for the parameters of the fitting line:
% we defined the search range in a way which enable detection of a
% replay which covers the whole track positions within 3 time bins
% wherever during the event; therefore the search space will be different
% for each event


if nargin > 0
    
    useReg = varargin{1};
    
end


if useReg
    
    [~, peakPos] = max(postPr, [], 1);
    peakPos      = peakPos';
    
    bins         = [ones(noTimeBins, 1) (1:noTimeBins)'];
    
    % Our regression can be based on only the peaks higher than an
    % arbitrary threshold
    
    sumAroundPeaks = zeros(noTimeBins, 1);
    for Tbin = 1:noTimeBins
        
        sumAroundPeaks(Tbin) = sum(postPr((max(1, peakPos(Tbin)-5):min(noPosBins, peakPos(Tbin)+5)), Tbin)); 
        
    end
    
    inclusionInx = sumAroundPeaks > 10/noPosBins;
    
    peakPos = peakPos(inclusionInx);
    bins    = bins(inclusionInx, :); 
    
    
    B            = bins\peakPos; % linear regression line (intercept and slope)
    
    initialRho   = B(1);
    initialSlope = B(2);
    
%     rhos   = initialRho - ceil(noPosBins/2) : 5: initialRho + ceil(noPosBins/2);
    rhos   = initialRho- noPosBins : 5: initialRho + noPosBins;
    slopes = initialSlope - noPosBins/5 : 2: initialSlope + noPosBins/5;
    
    
else
    
    rhos   = -(ceil(noTimeBins/3)-1) * noPosBins :10: ceil(noTimeBins/3) * noPosBins; 
    slopes = -ceil(noPosBins/2) :2: ceil(noPosBins/2); 
    
end



%%% calculating goodness of fit for each point in the seach space
% 
% neighborsSummed = 7; 
% postPr_summed = zeros(size(postPr));
% for ii = 1: noTimeBins
%     postPr_summed(:, ii) = conv(postPr(:, ii), ones(2*neighborsSummed+1, 1), 'same');
% end
% 


goodnessofFit = zeros(length(rhos), length(slopes));
for i = 1 : length(rhos)
     rho = rhos(i);
     for j = 1 : length(slopes)
         slope = slopes(j);
        
         goodnessofFit(i,j) = sumalongline(postPr, rho, slope); 
%          goodnessofFit(i,j) = sumalongline3(postPr_summed, rho, slope); 
         
     end
end


% 
% 
% lineParameters = [repelem(rhos', length(slopes), 1) repmat(slopes', length(rhos), 1)];
% lineParameters = mat2cell(lineParameters, ones(length(lineParameters), 1), 2);
% tic
% goodnessofFit = cellfun(@(x) sumalongline2(postPr, x), lineParameters);
% time2 = toc;




%%% finding the slope and rho (y intercept) resulting in the highest
%%% goodness of fit

maxgoodnessofFit = max(goodnessofFit(:));

[bestrhoInd, bestSlopeInd] = find(goodnessofFit == maxgoodnessofFit, 1, 'first');


bestrho = rhos(bestrhoInd);
bestSlope = slopes(bestSlopeInd);


%%% start and end points of the best fitting line

yBeg = bestrho + bestSlope*1; % first time bin
% yEnd = bestrho + bestSlope*(noTimeBins-1);
yEnd = bestrho + bestSlope*noTimeBins;

begPnt = [1 yBeg];
endPnt = [noTimeBins yEnd];

end

