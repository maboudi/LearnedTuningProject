function  sumProb = sumalongline3(postPrSummed, rho, slope)


noTimebins = size(postPrSummed,2);
noPositionBins = size(postPrSummed,1);

sumProb = 0;
% 
% middlePoints = ceil(rho + slope * (1:noTimebins)); 



for timebin = 1 : noTimebins
    
%     middlePoint = ceil(rho + slope * (timebin-1)); %% The position on the line in the current time bin
                                            %%% The sum of prhobability over
                                            %%% a vicinity of 30 cm (14 *
                                            %%% 2.5) arhound the middle point is calculated in the current time bin 
    
    middlePoint = ceil(rho + slope * timebin);
    
    if middlePoint > 0 && middlePoint < noPositionBins
        sumProb = sumProb + postPrSummed(middlePoint, timebin);
    else
        sumProb = sumProb + median(postPrSummed(:, timebin));
    end
    
    
end





end