function subEvent = timebinshuffle(actualEvents, eventnumber)


%   This function is intended to substitute bins within each event with a
%   random set of bins from the pool of events


    if size(actualEvents, 1) > size(actualEvents, 2)
        actualEvents = actualEvent';
    end

    cnctEvents = cell2mat(actualEvents);
    poolednoBins = size(cnctEvents, 2);

    
    noBins = size(actualEvents{eventnumber}, 2);
    rndgen = randi(poolednoBins, 1, noBins);

    subEvent = cnctEvents(:, rndgen);
         
end

