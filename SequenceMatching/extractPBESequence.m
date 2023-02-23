function [PBEsequence, smoothedPBE_sorted] = extractPBESequence(PBE)

% This function is intended to extract the sequence of activation of units
% within a PBE. This is done by smoothing the spike trains with a Gaussian
% kernel and sorting the units based on the peaks of smoothed curves.

% PBE: population burst event with 1 ms time bins


sigma = 20; % in ms
halfwidth = 2*sigma; 
gaussKernel = gausswindow(sigma, halfwidth);

activeUnits = find(sum(PBE, 2)); % the units active in the current population burst event
PBE_activeUnits = PBE(activeUnits, :);


nActUnits = numel(activeUnits);
smoothedPBE = zeros(nActUnits, size(PBE, 2));


for unit = 1:nActUnits; smoothedPBE(unit, :) = conv(PBE_activeUnits(unit, :), gaussKernel, 'same'); end

[~, peakTimeBins] = max(smoothedPBE, [], 2);


% PBE sequence

[~, sortInd] = sort(peakTimeBins, 'ascend');
PBEsequence = activeUnits(sortInd);

smoothedPBE_sorted = smoothedPBE(sortInd, :);
smoothedPBE_sorted = smoothedPBE_sorted ./repmat(max(smoothedPBE_sorted, [], 2), [1, size(smoothedPBE_sorted, 2)]);


end