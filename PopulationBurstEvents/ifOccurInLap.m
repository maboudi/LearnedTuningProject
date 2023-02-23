function PBElapIdx = ifOccurInLap(PBEs, laps)

nPBEs = size(PBEs, 1);
PBElapIdx = zeros(nPBEs, 1);

for ii = 1:nPBEs
    
    if isempty(find(isbetween(PBEs(ii, 1:3), laps), 1))
       PBElapIdx(ii) = 1;
    end

end

end


function isInsideInterval = isbetween(PBE, laps)

nLaps = length(laps);
isInsideInterval = zeros(nLaps, 1);

for ii = 1:nLaps
   
    if PBE(3) >= laps(ii, 1) & PBE(3) <= laps(ii, 2) 
        isInsideInterval(ii) = 1;
    end
    
end

end