function timeSwapEventsBinnedfiring = timeswap(eventsBinnedfiring, binDur)

% This function is intended to randomly permute the time bins

noEvents = size(eventsBinnedfiring, 1);

timeSwapEventsBinnedfiring = cell(noEvents, 1);

for evt = 1 : noEvents
    
    noBins = size(eventsBinnedfiring{evt, 2},2); % number of 20ms time bins
    
    rndgen = randperm(noBins); % random permutation indices for 20ms bins
    
    binSize = binDur*1000;

    indices = reshape(1:noBins*binSize, [binSize noBins]);
    rndgen2 = reshape(indices(:, rndgen), [noBins*binSize 1])'; % generate corresponding permutation indices for 1ms bins

    try
        timeSwapEventsBinnedfiring{evt, 1} = eventsBinnedfiring{evt, 1}(:, rndgen2);
    catch
        timeSwapEventsBinnedfiring{evt, 1} = [];
    end
    
    timeSwapEventsBinnedfiring{evt, 2} = eventsBinnedfiring{evt, 2}(:, rndgen);

    
end


end
