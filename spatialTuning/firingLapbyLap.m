function firingLapbyLap(spikeStruct, spatialTuningsA, spatialTuningsB, spatialInfoA, spatialInfoB, conslapsRatioA, conslapsRatioB, behavior, turningPeriods, fileinfo, directory)


spatialTuningsA = spatialTuningsA.smoothed;
spatialTuningsB = spatialTuningsB.smoothed;


tbegin = double(behavior.time(2,1));
Fs = fileinfo.Fs;
speedThresh = 10;
qclus = [1 2 3];

uniqUnits = unique(spikeStruct.unit);

bins = linspace(min(fileinfo.xyt2(:, 1)), max(fileinfo.xyt2(:, 1)), 20);
binCenters = bins(1:end -1)+ diff(bins(1:2))/2;

posbinSize = 2; % cm

for ii = 1:length(uniqUnits)
    
    currUnit = uniqUnits(ii);
    
    
    % direction A
    
    if ~isempty(spatialTuningsB)
        spikeInd = find(spikeStruct.unit == currUnit & spikeStruct.lap > 0 & ismember(spikeStruct.qclu, qclus) & ...
                        spikeStruct.speed > speedThresh  & mod(spikeStruct.lap, 2)== 0); % & spikeStruct.theta == 1
    else 
        spikeInd = find(spikeStruct.unit == currUnit & spikeStruct.lap > 0 & ismember(spikeStruct.qclu, qclus) & ...
                        spikeStruct.speed > speedThresh ); % & spikeStruct.theta == 1
    end
                
    turningSpikes = [];
    for jj = 1: size(turningPeriods, 1); turningSpikes = [turningSpikes; find(spikeStruct.t(spikeInd) >= turningPeriods(jj, 1) & spikeStruct.t(spikeInd) <= turningPeriods(jj, 2))]; end     
    
    turningSpikeInd = spikeInd(turningSpikes); 
    spikeInd(turningSpikes) = [];
    
                
    
    if isempty(ismember(unique(spikeStruct.qclu(spikeInd)), [1 2 3]))
       continue 
    end
    
    
    figure;
    set(gcf, 'unit', 'centimeter', 'position', [20 20 20 20])
   
    ax1 = subplot(7,2,[1 3 5 7 9]);
    plot(fileinfo.xyt2(:, 1), (fileinfo.xyt2(:, 2)-tbegin)/Fs, '.', 'color', [.5 .5 .5], 'markersize', 2); 
    hold on
    plot(spikeStruct.linearPos(spikeInd), (spikeStruct.t(spikeInd)-tbegin)/Fs, '.', 'color', [255 50 50]/255, 'markersize', 5)
    plot(spikeStruct.linearPos(turningSpikeInd), (spikeStruct.t(turningSpikeInd)-tbegin)/Fs, '.', 'color', 'g', 'markersize', 5)
    
    title({sprintf('unit %d', currUnit); 'direction A'},'fontsize', 10)
    
    ylabel('time(sec)')
    
    currYlim = ylim;
    currXlim = xlim;
    text(currXlim(1)+0.1*range(currXlim), currYlim(2)-0.06*range(currYlim), sprintf('spatial Info = %.2f\nconsistentlapsRatio = %.2f', spatialInfoA(currUnit), conslapsRatioA(currUnit)), 'color', 'k')

    ax2 = subplot(7,2,11);
    spkCnts = histc(spikeStruct.linearPos(spikeInd), bins); spkCnts(end) = [];
    bar(binCenters, spkCnts, 'edgecolor', 'none', 'facecolor', [255 50 50]/255) %,'facealpha', 0.5

    ylabel('total spike count')

    
    ax3 = subplot(7,2,13);
    
    linearPoscenters = (1:size(spatialTuningsA, 2))*posbinSize;
    
    fill([linearPoscenters fliplr(linearPoscenters)], [spatialTuningsA(currUnit, :) fliplr(zeros(size(spatialTuningsA(currUnit, :))))], [255 50 50]/255,'LineStyle','none')
    hold on
    plot(linearPoscenters, spatialTuningsA(currUnit, :),'color', [255 50 50]/255,'linewidth', 1);
%     alpha(0.5)
    
    ylabel('firing rate')
    xlabel('position (cm)')
    
    
    if ~isempty(spatialTuningsB)
    % direction B
    
    spikeInd = find(spikeStruct.unit == currUnit & spikeStruct.lap > 0 & ismember(spikeStruct.qclu, qclus) & ...
                    spikeStruct.speed > speedThresh  & mod(spikeStruct.lap, 2)== 1); % & spikeStruct.theta == 1
                
               
    turningSpikes = [];
    for jj = 1: size(turningPeriods, 1); turningSpikes = [turningSpikes; find(spikeStruct.t(spikeInd) >= turningPeriods(jj, 1) & spikeStruct.t(spikeInd) <= turningPeriods(jj, 2))]; end     
    
    turningSpikeInd = spikeInd(turningSpikes); 
    spikeInd(turningSpikes) = [];
    
    
    
    
    ax4= subplot(7,2,[2 4 6 8 10]);
    plot(fileinfo.xyt2(:, 1), (fileinfo.xyt2(:, 2)-tbegin)/Fs, '.', 'color', [.5 .5 .5], 'markersize', 2); 
    hold on
    plot(spikeStruct.linearPos(spikeInd), (spikeStruct.t(spikeInd)-tbegin)/Fs, '.', 'color', [51 150 255]/255, 'markersize', 5)
    plot(spikeStruct.linearPos(turningSpikeInd), (spikeStruct.t(turningSpikeInd)-tbegin)/Fs, '.', 'color', 'g', 'markersize', 5)

    title('direction B',  'fontsize', 10)
    
    currYlim = ylim;
    currXlim = xlim;
    text(currXlim(1)+0.1*range(currXlim), currYlim(2)-0.06*range(currYlim), sprintf('spatial Info = %.2f\nconsistentlapsRatio = %.2f', spatialInfoB(currUnit), conslapsRatioB(currUnit)), 'color', 'k')

    
%     set(ax4, 'yticklabel', [])

    ax5 = subplot(7,2,12);
    spkCnts = histc(spikeStruct.linearPos(spikeInd), bins); spkCnts(end) = [];
    bar(binCenters, spkCnts, 'edgecolor', 'none', 'facecolor', [51 150 255]/255) %,'facealpha', 0.5
    
%     set(ax5, 'yticklabel', [])
    
    ax6 = subplot(7,2,14);
    linearPoscenters = (1:size(spatialTuningsB, 2))*posbinSize;
    
    fill([linearPoscenters fliplr(linearPoscenters)], [spatialTuningsB(currUnit, :) fliplr(zeros(size(spatialTuningsB(currUnit, :))))], [51 150 255]/255,'LineStyle','none')
    hold on
    plot(linearPoscenters, spatialTuningsB(currUnit, :),'color', [51 150 255]/255,'linewidth', 1);
%     alpha(0.5)
    
%     ylabel('firing rate')
    xlabel('position (cm)')
%     set(ax6, 'yticklabel', [])

    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')
    linkaxes([ax1, ax4], 'y')
    
    else
        linkaxes([ax1 ax2 ax3], 'x')
    end
    
    
%     linkaxes([ax2, ax5], 'y')
%     linkaxes([ax3, ax6], 'y')
    

    
    filename = fullfile(directory, sprintf('LapbyLapfirings_unit%.3d', currUnit));
    print(gcf, filename,'-dpdf') 
        
    close all
end



