function [scoreZ, mcPval] = rawScoreCompare2null(rawScore, nullScoreDist)



med = median(nullScoreDist);
percentile95 = prctile(nullScoreDist, 95);

if rawScore ~= med
    
    if percentile95 == med
        
        %%% just assign some big value which its sign indicate whether the
        %%%raw score was to the right or to the left of the null dist(or value)
        scoreZ = 1e3 * sign(rawScore - med); 
        
    else
        scoreZ = (rawScore - med)/(percentile95 - med);
    end
    
else
    scoreZ = 0;
end

mcPval = (length(find(nullScoreDist >= rawScore)) + 1)/(length(nullScoreDist) + 1);

    
end