function idShuffleBinnedfiring = runUnitIDshuffle(spike, thetaBouts, qclus, activeUnits, binDur, Fs)


runSpikes = find(spike.theta == 1 & spike.lap > 0 & ismember(spike.qclu, qclus));
runSpikets = spike.t(runSpikes);
runUnits = spike.unit(runSpikes);
runqclus = spike.qclu(runSpikes);


tempSpike = struct('t', [], 'unit', []);
tempSpike.t = runSpikets; 

% reassiging new units IDs to the spikes

unique_units = unique(runUnits);

randgen = randperm(length(unique_units));
reassignedIDs = unique_units(randgen);


tempSpike.unit = zeros(length(runUnits), 1);
tempSpike.qclu = zeros(length(runqclus), 1);

for ii = 1 : length(unique_units)
    
    unitspikeidx = find(runUnits == unique_units(ii)); %% find all the spikes with an specific original unit ID
    
    tempSpike.unit(unitspikeidx) = reassignedIDs(ii); %% reassign the new random ID to the set of spikes
    tempSpike.qclu(unitspikeidx) = unique(runqclus(find(runUnits == reassignedIDs(ii))))';
    
end


tempDataBinned = timeBinning(thetaBouts, tempSpike, qclus, binDur, Fs);
tempDataBinned = justActiveUnits(tempDataBinned, activeUnits);

% remove the zero bins, tolerated zero bins is set to zero, minimum number
% of time bins is 4, at least four neurons should  fire in a theta event

idShuffleBinnedfiring = removeZeroBins(tempDataBinned, thetaBouts(:, 1), 0, 4, [], [], binDur, Fs); 
  
end
