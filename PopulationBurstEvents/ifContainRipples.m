function PBErippleIdx = ifContainRipples(PBEs, ripplePeriods)

nPBEs = size(PBEs, 1);

try
    peakRipples = ripplePeriods(:, 3);
catch
    peakRipples = mean(ripplePeriods, 2);
end

PBErippleIdx = zeros(nPBEs, 1);

for ii = 1:nPBEs
    
    if ~isempty(find(isbetween(peakRipples, PBEs(ii, 1:2)), 1))
        PBErippleIdx(ii) = 1;
    end

end

end


function isInsideInterval = isbetween(timePnts, interval)

ntps = length(timePnts);
isInsideInterval = zeros(ntps, 1);

for ii = 1:ntps
   
    if timePnts(ii) >= interval(1) && timePnts(ii) <= interval(2) 
        isInsideInterval(ii) = 1;
    end
    
end

end