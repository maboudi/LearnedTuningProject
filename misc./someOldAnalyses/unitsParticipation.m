function [participation, participantCells, partPCs, partnonPCs] = unitsParticipation(eventsBinnedfiring, activeUnits, RLtemplate, LRtemplate, fileinfo, longORshort, binDur)


noEvents = size(eventsBinnedfiring, 1);

%%% uni- and bidirectional units
noUnits = size(eventsBinnedfiring{1,1}, 1); 

allunits = 1 : noUnits;

leftonlyUnits = find(ismember(allunits, RLtemplate) & ~ismember(allunits, LRtemplate));
rightonlyUnits = find(ismember(allunits, LRtemplate) & ~ismember(allunits, RLtemplate));
bidirectionalUnits = find(ismember(allunits, RLtemplate) & ismember(allunits, LRtemplate));
notplaceUnits = find(~ismember(allunits, RLtemplate) & ~ismember(allunits, LRtemplate) & ismember(allunits, activeUnits));



%%% firing status of the units across population burst events

concatEvents = cell2mat(eventsBinnedfiring(:,2)');
meanFiring = mean(concatEvents,2);
meanFiring = meanFiring/binDur;
firing95prctile = prctile(meanFiring, 95);
 

figure; %('Visible','off')

fHand = subplot(2,2,1);

hold(fHand, 'on')

for unit = 1 : noUnits
    
    if ismember(unit, leftonlyUnits)
        singleBarColor = 'r';
    elseif ismember(unit, rightonlyUnits)
        singleBarColor = 'b';
    elseif ismember(unit, bidirectionalUnits)
        singleBarColor = 'm';
    elseif ismember(unit, notplaceUnits)
        singleBarColor = 'k';
    else
        singleBarColor = 'g';
    end
    
    bar(unit, meanFiring(unit), 'parent', fHand, 'facecolor', singleBarColor, 'edgecolor','none')
end

line([1 noUnits], [firing95prctile firing95prctile], 'linewidth', 2)

hold off

xlabel('Unit ID', 'fontsize', 20)
ylabel('Average Firing Rate(Hz)', 'fontsize', 20)
set(gca,'fontsize', 16);





%%%% calculating the number of events that each unit participated;
%%% And also the number of place cells and non-place cells (only putative
%%%% pyramiodal) which has participated in each event


participation = zeros(noUnits, noEvents); %%% participation status of each unit in each event (noUnits * noEvents)

participantCells = cell(1, noEvents);
noPartCells = zeros(noEvents, 1);

partPCs = cell(1, noEvents);
partnonPCs = cell(1, noEvents);

noPartPCs = zeros(noEvents, 1);
noPartnonPCs = zeros(noEvents, 1);



for evt = 1 : noEvents
    
    participation(:,evt) = sum(eventsBinnedfiring{evt,2},2) > 0;
    
    %%% all participant units in the current event
    
    participantCells{evt} = find(participation(:,evt));
    
    noPartCells(evt) = length(participantCells{evt}); %% the total number
    
    
    %%% all place units and non-place units participating in the current
    %%% events, separately
    
    partPCs{evt} = participantCells{evt}(ismember(participantCells{evt}, unique([RLtemplate; LRtemplate]))); %% participant pc cells
    partnonPCs{evt} = participantCells{evt}(~ismember(participantCells{evt}, unique([RLtemplate; LRtemplate])));
    
    noPartPCs(evt) = length(partPCs{evt});
    noPartnonPCs(evt) = length(partnonPCs{evt});
    
end

participationRatio = mean(participation, 2) * 100; %%% participation ratio of each unit
PR95prctile = prctile(participationRatio, 95); %% paricipation ratio (PR) 95 percentile

%%%%
fHand = subplot(2,2,3);

hold(fHand, 'on')


for unit = 1 : noUnits
    
    if ismember(unit, leftonlyUnits)
        singleBarColor = 'r';
    elseif ismember(unit, rightonlyUnits)
        singleBarColor = 'b';
    elseif ismember(unit, bidirectionalUnits)
        singleBarColor = 'm';
    elseif ismember(unit, notplaceUnits)
        singleBarColor = 'k';
    else
	singleBarColor = 'g';
    end
    
    bar(unit, participationRatio(unit), 'parent', fHand, 'facecolor', singleBarColor, 'edgecolor','none')
end

line([1 noUnits], [PR95prctile PR95prctile], 'linewidth', 2)

hold off

xlabel('Unit ID', 'fontsize', 20)
ylabel('Event Participation(%)', 'fontsize', 20)
set(gca,'fontsize', 16);



subplot(2,2,2)

[h1, bins] = hist(noPartPCs, 50);
h2 = hist(noPartnonPCs, bins);

hold on
bar(bins, h1, 'facecolor', 'm', 'edgecolor', 'none')
% hold on
bar(bins, h2, 'facecolor', 'none', 'edgecolor', 'k');

hold off

xlabel('Number of participant cells', 'fontsize', 20)
ylabel('Number of events', 'fontsize', 20)
set(gca,'fontsize', 16);


subplot(2,2,4)

noPartPCsRatio = (noPartPCs ./ noPartCells) * 100;  

[h1, bin] = hist(noPartPCsRatio, 50);
bar(bin, h1, 'facecolor', 'm', 'edgecolor', 'none')

xlabel('Ratio of participant Place cells', 'fontsize', 20)
ylabel('Number of events', 'fontsize', 20)
set(gca,'fontsize', 16);

% currDir = pwd;
% Folderbase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/Population Burst Events'];

% saveas(gcf, [Folderbase '/' 'UnitsEventParticipation.fig'])
% save([Folderbase '/ParticipationInfo.mat'], 'participation', 'participantCells', 'partPCs', 'partnonPCs')

end
