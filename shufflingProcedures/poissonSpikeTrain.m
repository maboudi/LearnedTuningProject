function    poissonEventsBinnedfiring = poissonSpikeTrain(eventsBinnedfiring, binDur)

% This function is intended to generate binned spike count for each unit. 
% The mean firing of each unit is calculated over the pool of events, which gives us the lambda. 
% The lambda is used to generate new spike trains



binSize = binDur * 1000;

noEvents = size(eventsBinnedfiring, 1);


%%% calculating the mean firing of the neurons across all the events (we
%%% will need this to generate Poisson simulated data

temp = cell2mat(eventsBinnedfiring(:,2)'); %% concatenate all the events
meanFiring = mean(temp,2);
meanFiring = meanFiring / binDur;

dt = 1/1000; %s


%%% now we define a poisson simulated data corresponding to each event
poissonEventsBinnedfiring = cell(noEvents, 2);

for evt = 1 : noEvents
    
    noUnits = size(eventsBinnedfiring{evt, 2}, 1); 
    noBins = size(eventsBinnedfiring{evt, 2}, 2);
    
    noMiliBins = noBins * binSize; %% number of milisecond bins
    
    edges = 0 : binSize : noBins*binSize;
    
%     poissonEventsBinnedfiring{evt, 1} = rand(noUnits, noMiliBins) < repmat((meanFiring  * dt), [1 noMiliBins]);
    
    
    for unit = 1 : noUnits
        
        poissonEventsBinnedfiring{evt, 1}(unit, :) = rand(1, noMiliBins) < meanFiring(unit) * dt;
        
        temp = histc(find(poissonEventsBinnedfiring{evt, 1}(unit, :)) , edges);
        poissonEventsBinnedfiring{evt, 2}(unit, :) = temp(:, 1:end-1);
    end
   
end

 
end 