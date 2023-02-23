function [postPr, fitness, bestSlope, bestRo, startPnt, endPnt] = replayevaluation(data)

     postPr = baysDecoder(data, LambdaMat, 0.015);
     
     noPosBins = size(LambdaMat, 2);
     noTimeBins = size(data, 2);
     
     ros = [-noPosBins : 2*noPosBins];
     slopes = -noPosBins : noPosBins;
     
     goodnessfit = zeros(length(ros), length(slopes));
     for i = 1 : length(ros)
         ro = ros(i);
         for j = 1 : length(slopes)
             slope = slopes(j);
             goodnessfit(i,j) = fittingline(postPr, ro, slope);
         end
     end
     
     [fitness roSlope] = max(goodnessfit(:));
     bestSlope = slopes(ceil(roSlope/length(ros)));
     if ceil(roSlope/length(ros)) == roSlope/length(ros)
        bestRo = ros(end);
     else
        bestRo = ros(roSlope - floor(roSlope/length(ros))*length(ros));
     end

     ySt = bestRo;
     yEnd = bestRo + bestSlope*(noTimeBins-1);
     
     startPnt = [1 ySt];
     endPnt = [noTimeBins yEnd];
     
end

