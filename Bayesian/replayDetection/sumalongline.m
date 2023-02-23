function sumProb = sumalongline(postPr, rho, slope)

noTimebins = size(postPr,2);
noPositionBins = size(postPr,1);

sumProb = 0;


for timebin = 1 : noTimebins
    
%     middlePoint = ceil(rho + slope * (timebin-1)); %% The position on the line in the current time bin
                                            %%% The sum of prhobability over
                                            %%% a vicinity of 30 cm (14 *
                                            %%% 2.5) arhound the middle point is calculated in the current time bin 
    
    middlePoint = ceil(rho + slope * timebin);
                                            
                                            
    if  middlePoint+7 >= 1 && middlePoint-7 <= noPositionBins %%% considering only the timepoints with a vicinity of middlepoint overlapping the track poisitions
        
        sumProb = sumProb + sum(postPr(max(1,(middlePoint-7)):min((middlePoint+7), noPositionBins), timebin));

    else
        continue
    end  
end

% timebins = 1:noTimebins;
% middlePoints = ceil(rho + slope * timebins);
% 
% sumOverArray = repmat(middlePoints, [15, 1]) - repmat(transpose(-7:7), [1, noTimebins]);
% 
% sumOverArray(sumOverArray < 1) = nan;
% sumOverArray(sumOverArray > noPositionBins) = nan;
% 
% postPr_1D = postPr(:);
% 
% sumOverArray_1D = sumOverArray + repmat([0:noTimebins-1]*noPositionBins, [15, 1]);
% sumOverArray_1D = sumOverArray_1D(:);
% 
% sumOverArray_1D(isnan(sumOverArray_1D)) = [];
% 
% sumProb = sum(postPr_1D(sumOverArray_1D));



noFiringTimeBins = length(find(sum(postPr, 1) > 0.9)); % in the non-firing time bins the sum is almost zero

sumProb = sumProb/noFiringTimeBins;


end