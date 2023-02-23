function [eventsBinnedfiring2, newEventbeg, newEventend, originalIndex]= removeZeroBins(eventsBinnedfiring, eventBeg, permittedNumofZeroBins, minNoBinsperSegment, fileinfo, whichPart, binDur, Fs)


NumofEvents = size(eventsBinnedfiring, 1); %% original events
binSize = double(binDur * 1000); 



newEvtInd = 0; %% index of the new event (which is resulted by segmenting the orignial events)


presetNumofNewEvents = NumofEvents * 4; %% we preset the number of new events as the number of original events time 4

originalIndex = zeros(presetNumofNewEvents , 2);
eventsBinnedfiring2 = cell(presetNumofNewEvents, 2);
newEventbeg = zeros(presetNumofNewEvents, 1);
newEventend = zeros(presetNumofNewEvents, 1);

for evt = 1 : NumofEvents
    
    eventData = eventsBinnedfiring{evt, 2};
    
    binSumFiring = sum(eventData, 1); %% sum of the firing over the units within each bin
    
    if isempty(find(binSumFiring, 1)) %% if no there is no firing bins then eliminate that event
        continue
    end
    
    
    %%% finding the starts and the ends of the firing segments
    
    firingBins = find(binSumFiring); 
    
    newEvtbeg = firingBins(find([0 diff(firingBins)] > 1));
    newEvtbeg = [firingBins(1) newEvtbeg];
    
    newEvtend = firingBins(find([0 diff(firingBins)] > 1) - 1);
    newEvtend = [newEvtend firingBins(end)];
    
    
    %%% merge the successive events in case there is just one zero bin
    %%% separating them

    for seg = 2 :length(newEvtbeg)
        
        if newEvtbeg(seg) - newEvtend(seg-1) - 1 <= permittedNumofZeroBins
           newEvtbeg(seg) = 0;
           newEvtend(seg-1) = 0;
        end
        
    end
    
    newEvtbeg = newEvtbeg(find(newEvtbeg));
    newEvtend = newEvtend(find(newEvtend));

    segNoBins = newEvtend - newEvtbeg + 1;
    
    
    for seg = 1 : length(newEvtbeg)
        
        currData = eventData(:, newEvtbeg(seg):newEvtend(seg));
        NumFiringUnits = length(find(sum(currData, 2)));
        
        if segNoBins(seg) >= minNoBinsperSegment && NumFiringUnits >= 5
            
           newEvtInd = newEvtInd + 1;
           originalIndex(newEvtInd,:) = [newEvtInd evt]; %%% save index of the new event with the index of original event
           
        else
            continue
        
        end
        
        eventsBinnedfiring2{newEvtInd, 2} = eventData(:, newEvtbeg(seg):newEvtend(seg));
        
        try
        eventsBinnedfiring2{newEvtInd, 1} = eventsBinnedfiring{evt, 1}(:, (newEvtbeg(seg)-1)*binSize+1 : newEvtend(seg)*binSize);
        catch
            (newEvtbeg(seg)-1)*binSize+1
            newEvtend(seg)*binSize
        end
        
        newEventbeg(newEvtInd) = eventBeg(evt) + (newEvtbeg(seg)-1)*binDur*Fs+1;
        newEventend(newEvtInd) = eventBeg(evt) + newEvtend(seg)*binDur*Fs;
    end
end

numofNewEvents = newEvtInd;


%%% removing the zero (empty) elements (cells) from the ends of the arrays
%%% with unknown size at initialization 

eventsBinnedfiring2 = eventsBinnedfiring2(1:numofNewEvents, :);
newEventbeg = newEventbeg(1:numofNewEvents);
newEventend = newEventend(1:numofNewEvents);
originalIndex = originalIndex(1:numofNewEvents, :);


if ~isempty(fileinfo)
    currDir = pwd;
    Folderbase = [currDir '/' fileinfo.name '/data part' num2str(whichPart) '/Population Burst Events']; %%'/data part' num2str(whichPart)
    save([Folderbase '/' 'ZeroRemovedEvents.mat'], 'originalIndex', 'newEventbeg', 'newEventend')

    %%% storing the new events

    Filename = [Folderbase '/' fileinfo.name '.fbe.evt'];
    MakeEvtFile([newEventbeg newEventend], Filename, {'beg', 'end'}, Fs, 1)
end


end