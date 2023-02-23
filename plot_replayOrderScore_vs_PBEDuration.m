


PBEInfo_PRE  = PBEInfo_selectedEvents(strcmp({PBEInfo_selectedEvents.epoch}, 'pre'));
PBEInfo_POST = PBEInfo_selectedEvents(strcmp({PBEInfo_selectedEvents.epoch}, 'post'));
PBEInfo_RUN  = PBEInfo_selectedEvents(strcmp({PBEInfo_selectedEvents.epoch}, 'run'));




calReplayScore = 0;
PBEInfo_PRE = BDreplayDetect_separateDirs(PBEInfo_PRE, spikes, fileInfo, binDur, calReplayScore);

PBEInfo_RUN = BDreplayDetect_separateDirs(PBEInfo_RUN, spikes, fileInfo, binDur, calReplayScore);

PBEInfo_POST = BDreplayDetect_separateDirs(PBEInfo_POST, spikes, fileInfo, binDur, calReplayScore);




%
temp.PRE = PBEInfo_PRE;
temp.RUN = PBEInfo_RUN;
temp.POST = PBEInfo_POST;

epochNames = {'PRE'; 'RUN'; 'POST'};


figure;
for epoch = 1:3
    
    currEpoch = epochNames{epoch};
    
    radonMixDir.(currEpoch)  = zeros(numel(temp.(currEpoch)), 1);
    weightedCorr.(currEpoch) = zeros(numel(temp.(currEpoch)), 1);
    
    for pbe = 1:numel(temp.(currEpoch))
        radonMixDir.(currEpoch)(pbe)  = temp.(currEpoch)(pbe).radonIntegral.mix;
        weightedCorr.(currEpoch)(pbe) = temp.(currEpoch)(pbe).weightedCorr.mix;
    end
    
    replayDurations.(currEpoch)   = [temp.(currEpoch).endT] - [temp.(currEpoch).startT];
    replayOrderScores.(currEpoch) = [temp.(currEpoch).replayOrderScore];
    replayOrderScores_prctile.(currEpoch) = [temp.(currEpoch).replayOrderScore_prctile];
    
    
    ax(epoch) = subplot(1,3,epoch);
    scatter(replayDurations.(currEpoch), replayOrderScores.(currEpoch), 5, replayOrderScores_prctile.(currEpoch), 'filled')
    colormap('jet')
    xlabel('Replay duration (s)')
    ylabel('Replay order')
    
    
    if epoch == 3
        ch = colorbar;
        colorTitleHandle = get(ch,'Title');
        set(colorTitleHandle ,'String', 'percentile');
    end
    
end
linkaxes(ax, 'y')

suptitle(fileInfo.sessionName)




% figure; 
% 
% ax(1)=subplot(3,1,1);
% hist([PBEInfo_PRE.replayOrderScore], 30)
% title('PRE')
% ylabel('number of PBEs')
% 
% ax(2)= subplot(3,1,2);
% hist([PBEInfo_RUN.replayOrderScore], 30)
% title('RUN')
% ylabel('number of PBEs')
% 
% ax(3) = subplot(3,1,3);
% hist([PBEInfo_POST.replayOrderScore], 30)
% title('POST')
% xlabel('replay order score')
% ylabel('number of PBEs')
% 
% linkaxes(ax, 'x')


