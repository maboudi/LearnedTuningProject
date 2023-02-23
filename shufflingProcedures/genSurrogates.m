
function BinnedEvts = genSurrogates(BinnedEvts)


%% Generating shuffle surrogate PBE datasets

% generate poisson simulated dataset

binDur = 0.02;
% noEvents = size(BinnedEvts, 1);

BinnedEvts.p = poissonSpikeTrain(BinnedEvts.data, binDur); 
% cmprPoisson2Raw(poissonBinnedEvts, BinnedEvts, longORshort, fileinfo, binDur); %%% compare the poisson simulated data with the actual data distribution

% generate time swap dataset (coherent shuffle within each event, keeping the cofirings the same)

BinnedEvts.ts = timeswap(BinnedEvts.data, binDur);
 
% time bins temporal shuffle (independently for each unit, incoherent shuffle)

% temporalBinnedEvts = cell(noEvents, 1);
% for evt = 1 : noEvents
%     
%     currEvt = BinnedEvts(evt, :);
%     
%     randShift = randi(size(currEvt{2}, 2)-1, 1, size(currEvt{2}, 1)); % generate shift amount for each row
%     
%     
%     % Because we are doing row cycle shuffle, we need to transpose the data
%     % matrix and after doing column cycle shuffles transpose it again
%     
%     temporalBinnedEvts{evt, 2} = column_cycle_shuffle(currEvt{2}', randShift)';
%     temporalBinnedEvts{evt, 1} = column_cycle_shuffle(currEvt{1}', randShift*20)'; % note the shift here
%     
% end

BinnedEvts.pts       = cell(size(BinnedEvts.data));
BinnedEvts.pts(:, 2) = genTimeSwap(BinnedEvts.data(:, 2));


end
