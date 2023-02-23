function jitterBinnedfiring = runSpikeJitter(spike, thetaBouts, jitterWinDur, qclus, activeUnits, binDur, Fs)


runSpikes = find(spike.theta == 1 & spike.lap > 0 & ismember(spike.qclu, qclus));
runSpikets = spike.t(runSpikes);
runUnits = spike.unit(runSpikes);
runqclus = spike.qclu(runSpikes);

numofSpikes = length(runSpikets);

tempSpike = struct('t', [], 'unit', [], 'qclu', []);
tempSpike.unit = runUnits;
tempSpike.qclu = runqclus;

tempSpike.t = runSpikets + sign(rand(numofSpikes, 1) - 0.5) .* rand(numofSpikes, 1) .* jitterWinDur; 


tempDataBinned = timeBinning(thetaBouts, tempSpike, qclus, binDur, Fs);
tempDataBinned = justActiveUnits(tempDataBinned, activeUnits);

% remove the zero bins, tolerated zero bins is set to zero, and minimum number
% of time bins is 4

jitterBinnedfiring = removeZeroBins(tempDataBinned, thetaBouts(:, 1), 0, 4, [], [], binDur, Fs); %% (:,:, jdata)
  
end
