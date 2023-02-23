function [varCollect, var_of_interest_mean, var_of_interest_sem, counts, binCenters] = timebinning_sleepEvolution(var_of_interest, peakT, binDur, stepDur, period)

totalT = period(2) - period(1);
nBins  = floor((totalT - binDur)/stepDur) + 1;

startT     = period(1) + (0:nBins-1)*stepDur;
endT       = startT + binDur;
binCenters = startT + binDur/2;

varCollect = cell(nBins, 1);
var_of_interest_mean = zeros(nBins, 1);
var_of_interest_sem = zeros(nBins, 1);
counts = zeros(nBins, 1);

for ii = 1:nBins
    
   idx = peakT > startT(ii) & peakT <= endT(ii);
   
   varCollect{ii} = var_of_interest(idx);
   var_of_interest_mean(ii) = nanmean(varCollect{ii});
   var_of_interest_sem(ii)  = nanstd(varCollect{ii})/numel(find(idx));
   counts(ii) = numel(varCollect{ii});
    
end

end