clear; clc;
load 2006-6-07_11-26-53.spikeII.mat
pos = load('2006-6-07_11-26-53.whl');
psf = 60; % position sampling frequency (Hz)
pix2cm = 8/2.78;
ConvPos = pos * pix2cm;

cpx = ConvPos(:,1);
cpy = ConvPos(:,2);

runWin = [ones(15,1)/15];
cpx = filtfilt(runWin,1,cpx);
cpy = filtfilt(runWin,1,cpy);

rcurr = sqrt(cpx.^2 + cpy.^2);

rprv = [nan; rcurr];
rnxt = rcurr(2:end);

tpnts = length(rcurr);
speed = zeros(tpnts,1);

speed(2:tpnts-1) = (rnxt(2:tpnts-1) - rprv(2:tpnts-1))/(2/psf); 
speed(1) = (rnxt(1) - rcurr(1))/(1/psf);
speed(tpnts) = (rcurr(tpnts) - rprv(tpnts))/(1/psf); %% 

speed = abs(speed);
stopPosInd = find(speed < 5);
bounds = find(diff(stopPosInd) > 1);
bndStr = [stopPosInd(1); stopPosInd(bounds+1)];
bndEnd = [stopPosInd(bounds); stopPosInd(end)];

speedfil = zeros(length(speed)*1667, 1); %% 0.01 ms binning (then we will  convert it to milisecond bins)
for i = 1 : length(bndStr)
    speedfil(((bndStr(i)-1)*1667 + 1) : bndEnd(i) * 1667) = 1;
end

speedfil = speedfil(1 : 100*floor(length(speedfil)/100));
temp = reshape(speedfil,[100, length(speedfil)/100]);
speedfil = mean(temp)';


Fs = 32552;

spikeTrain = spike.t(floor(1:length(spike.t)/1));
spikePos = interp1(0:1/60:(length(pos)-1)/60, cpx, spikeTrain/Fs);
spikeSpeed = interp1(0:1/60:(length(pos)-1)/60, speed, spikeTrain/Fs);

trainStr = min(spikeTrain);
trainEnd = max(spikeTrain);


spkTrnDur = trainEnd;

%%% time binnig of the spike train, bin size = 1 ms

binSize = 0.001*Fs;
noBins = floor(spkTrnDur / binSize);
        
binEdges = zeros(noBins + 1,1);
binEdges(1) = 1; %trainStr
for edge = 1 : noBins
    binEdges(edge+1) = edge * binSize;
end

[spikeBinned, binInd] = histc(spikeTrain, binEdges);


%%% spike density function

sigma = 20; %% in milisecond as we had use 1 milisecond time bins
halfwidth = 3 * sigma;

smoothwin = gausswindow(sigma, halfwidth);
temp = conv(spikeBinned, smoothwin);
spkDensity = temp(halfwidth+1 : end-halfwidth);

% figure;
% hist(spkDensity, 100)
stopBins = find(speedfil(1:length(spkDensity)) == 1);
rateSD = std(spkDensity(stopBins)); %%(stopBins)
rateMean = mean(spkDensity(stopBins)); %%(stopBins)

highRateBins = find(spkDensity > (rateMean + 2*rateSD));
boundaries = highRateBins(find(diff(highRateBins) > 1)+1);

epochStr(1) = highRateBins(1);
for epoch = 2 : length(boundaries+1)
    epochStr(epoch) = boundaries(epoch-1);
end

epochEnd(length(boundaries+1)) = highRateBins(end); %% the last bin of the last epoch
for epoch = 1 : length(boundaries) %% for the first epoch to the one before the end
    epochEnd(epoch) = boundaries(epoch)-1;
end

%%% find the bin with peak rate in each of the epoches
for epoch = 1 : length(epochStr)
    [maxpeak(epoch), Ind] = max(spkDensity(epochStr(epoch):epochEnd(epoch)));
    
    peakInd(epoch) = Ind + epochStr(epoch) - 1;
end

Evts = [epochStr; epochEnd; peakInd]';
eventSpeed = interp1(0:1/60:(length(pos)-1)/60, speed, peakInd/1e3);
stopEvents = find(eventSpeed < 5);



Evts = Evts(stopEvents, :);
evtHighestMUARate = (maxpeak(stopEvents) - rateMean)/rateSD;

%%% determine the firing rate rising and falling edges: the criterion used
%%% here is the mean firing rate across all the spike train

halfmaxlength = 500;

er = zeros(1, size(Evts, 1));
for epoch = 1 : size(Evts, 1)
    
    if Evts(epoch, 3) < halfmaxlength
        startcheckpoint = 1;
    else
        startcheckpoint = Evts(epoch, 3) - halfmaxlength;
    end
    
    if Evts(epoch, 3) > (noBins - halfmaxlength)
        endcheckpoint = noBins;
    else
        endcheckpoint = Evts(epoch, 3) + halfmaxlength;
    end
    
    meanCutpoint = find(spkDensity(startcheckpoint : Evts(epoch, 3)) < rateMean, 1, 'last');
    
    difftemp = [0; diff(spkDensity(startcheckpoint : Evts(epoch, 3)))];
    changesign = [0; difftemp(1:end-1).*difftemp(2:end)];
    zerodiff = find(changesign < 0 & spkDensity(startcheckpoint : Evts(epoch, 3)) < rateMean, 1, 'last')-1;
    
    if ~isempty(meanCutpoint) & ~isempty(zerodiff)
        bndPnt1 = min(meanCutpoint, zerodiff);
    elseif ~isempty(meanCutpoint)
        bndPnt1 = meanCutpoint;
    elseif ~isempty(zerodiff)
        bndPnt1 = zerodiff;
    else
        er(epoch) = 1;
    end
    
    if er(epoch) ~= 1
        
        risingEdge(epoch) = Evts(epoch, 3) - (halfmaxlength-bndPnt1)-1;
    else
        risingEdge(epoch) = Evts(epoch, 3) - halfmaxlength;
    end
    
    meanCutpoint = find(spkDensity(Evts(epoch, 3) : endcheckpoint) < rateMean, 1, 'first');
    
    difftemp = [0; diff(spkDensity(Evts(epoch, 3) : endcheckpoint))];
    changesign = [0; difftemp(1:end-1).* difftemp(2:end)];
    zerodiff = find(changesign < 0 & spkDensity(Evts(epoch, 3) : endcheckpoint) < rateMean, 1, 'first')-1;
    
    if ~isempty(meanCutpoint) & ~isempty(zerodiff)
        bndPnt2 = max(meanCutpoint, zerodiff);
    elseif ~isempty(meanCutpoint)
        bndPnt2 = meanCutpoint;
    elseif ~isempty(zerodiff)
        bndPnt2 = zerodiff;
    else
        er(epoch) = 1;
    end
    
    if er(epoch) ~= 1
        fallingEdge(epoch) = Evts(epoch, 3) + bndPnt2 - 1;
    else
        fallingEdge(epoch) = Evts(epoch, 3) + halfmaxlength;
    end
    
    epochLen(epoch) = fallingEdge(epoch) - risingEdge(epoch);
end

%%% finding epochs (events) which satisfy duration criterion
accepted_epochs = find(er ~= 1 & epochLen >= 60 & epochLen <= 500);
primaryEvts = [risingEdge(accepted_epochs) ; fallingEdge(accepted_epochs) ; Evts(accepted_epochs, 3)']';
evtHighestMUARate = evtHighestMUARate(accepted_epochs);

%%% detecting edges which are repeated and should be removed, because the
%%% edges are need to be non-decreasing

risefall_temp = [primaryEvts(:,1) primaryEvts(:,2)]';
edges = risefall_temp(:);

repeated_edges = find([edges(1); diff(edges)] < -30); %%%% confusing a bit, work more on it
secondaryEvts = primaryEvts;
secondaryEvts(floor(repeated_edges/2)+1 , :) = [];
evtHighestMUARate(floor(repeated_edges/2)+1) = [];

EvtsLength = secondaryEvts(:,2) - secondaryEvts(:,1);
%%% distribution of lengths of population burst events
figure;
hist(EvtsLength, 50)

eventPos = interp1(0:1/60:(length(pos)-1)/60, cpx, secondaryEvts(:,3)/1e3);


% 
% risefall_temp = [secondaryEvts(:,1) secondaryEvts(:,2)]';
% edges = risefall_temp(:);
% 
% spikehighpop = zeros(length(spike.t),1);
% [temp spikeInd] = histc(spike.t, edges*binSize);
% 
% for spk = 1 : length(spike.t)
%     if ceil(spikeInd(spk)/2) ~= spikeInd(spk)/2 
%         spikehighpop(spk) = 1;
%     end
% end



% rippleStr = spike.t(find([spikehighpop(1); diff(spikehighpop)] == 1));
% rippleEnd = spike.t(find([spikehighpop(1); diff(spikehighpop)] == -1)-1);

% secondaryEvtsInFrames = (secondaryEvts-1)*binSize+1;
% save('spikeEvent.mat', 'secondaryEvtsInFrames')
% 
% MakeEvtFile(secondaryEvtsInFrames, '2006-6-07_11-26-53.stp.evt', {'beg', 'end', 'peak'}, Fs,1)
% 


figure;

% subplot(5,1,[1 4])
% 
% hold on
% plot(spikeTrain/Fs , spike.cluster(floor(1:length(spike.t)/1)), '.', 'color', 'k', 'markersize', 5);


% for i = 1 : length(spikeTrain)
%     plot([spikeTrain(i) spikeTrain(i)]/Fs, [spike.cluster(i) - 0.2 spike.cluster(i) + 0.2], 'color','k')
% end

% yL = get(gca,'YLim');
% for epoch = 1:length(epochStr) 
% 
%     if ismember(epoch, slowepochInd) & er(epoch) ~= 1 & epochLen(epoch) >= 70
%         
%         line([risingEdge(epoch) risingEdge(epoch)]*binSize/Fs , yL, 'color', 'r', 'linestyle',':');
%         line([fallingEdge(epoch) fallingEdge(epoch)]*binSize/Fs , yL, 'color', [.7 .7 .7], 'linestyle',':');
%         
%         
%        p = patch([risingEdge(epoch) fallingEdge(epoch) fallingEdge(epoch) risingEdge(epoch)]*binSize/Fs,[0 0 yL(2) yL(2)], ...
%            'r','edgecolor', 'none');
%        set(p,'FaceAlpha', 0.1)
%     end
% end
% hold off
% xlim([2500000 3150000]/Fs)
% set(gca,'xtick', []);
% ylabel('cell ID', 'fontsize', 20)


% subplot(5,1,5)
hold on
t = (1 : binSize : 1 + noBins*binSize)/Fs;
%t = (trainStr : binSize : trainStr + noBins*binSize)/Fs;
plot(t, spkDensity, 'linewidth', 1, 'color','k')
% plot(t, rateMean*10, 'r', 'linewidth', 2)
yL = get(gca,'YLim');

for epoch = 1 : size(secondaryEvts,1)
    line([(secondaryEvts(epoch,1)-1)*binSize+1 (secondaryEvts(epoch,1)-1)*binSize+1]/Fs , yL, 'color', 'r', 'linestyle',':');
    line([(secondaryEvts(epoch,2)-1)*binSize (secondaryEvts(epoch,2)-1)*binSize]/Fs , yL, 'color', [.7 .7 .7], 'linestyle',':');


    p = patch([(secondaryEvts(epoch,1)-1)*binSize+1 (secondaryEvts(epoch,2)-1)*binSize (secondaryEvts(epoch,2)-1)*binSize (secondaryEvts(epoch,1)-1)*binSize+1]/Fs,[0 0 yL(2) yL(2)], ...
       'r','edgecolor', 'none');
    set(p,'FaceAlpha', 0.1)
end


%%% find the sharp wave/ ripple periods (is not precise)

rippleStr = spike.t(find([spike.ripple(1); diff(spike.ripple(floor(1:length(spike.t)/1)))] == 1));
rippleEnd = spike.t(find([spike.ripple(1); diff(spike.ripple(floor(1:length(spike.t)/1)))] == -1));
noRipple = length(rippleStr);
ripplePos = interp1(0:1/60:(length(pos)-1)/60, cpx, rippleStr/Fs);

% 
% overlapping_ripple = zeros(1, noRipple);
% for k = 1 : noRipple
%     p = patch([rippleStr(k) rippleEnd(k) rippleEnd(k) rippleStr(k)]/Fs,[0 0 yL(2) yL(2)],'b', 'edgecolor','none');
%     set(p, 'FaceAlpha', 0.1)
%     
%     ripple_checked = 0;
%     for epoch = accepted_epochs
%         if ((rippleStr(k) >= risingEdge(epoch)*binSize & rippleStr(k) <= fallingEdge(epoch)*binSize) | ...
%                 (rippleEnd(k) >= risingEdge(epoch)*binSize & rippleEnd(k) <= fallingEdge(epoch)*binSize)) & ripple_checked == 0
%             overlapping_ripple(k) = 1;
%             ripple_checked = 1;
%         end
%     end      
% end

%%% Number of epochs common between high polpulation rate epochs and SWR
%%% periods
% noOverlap = length(find(overlapping_ripple));

plot(t, ones(1,length(t))*rateMean, '--r', 'linewidth', 2)
plot(t, ones(1,length(t))*(rateMean + 2*rateSD), 'r', 'linewidth', 2)

hold off
% xlim([1000000 3150000]/Fs)
ylabel('firing rate (kHz)', 'fontsize', 20)
xlabel('time(sec)', 'fontsize', 20)



%%% ploting the raster of population firing related to one event for
%%% illustration 
% figure;
% hold on
% insetRasterInd = find(spike.t <= fallingEdge(29)*binSize & spike.t >= risingEdge(29)*binSize);
% insetRaster = spike.t(insetRasterInd);
% for ii = 1 : length(insetRaster)
%     plot([insetRaster(ii) insetRaster(ii)]/Fs * 1000, [spike.cluster(insetRasterInd(ii)) - 0.2 spike.cluster(insetRasterInd(ii)) + 0.2], 'color', 'k')
% end
% 
% hold off
% xlabel('time(msec)', 'fontsize', 16)

figure;
subplot(5,1,1:3)
plot(spikePos, spikeTrain/Fs, '.','color', 'k','markersize', 0.5)
hold on
plot(spikePos(spikeSpeed < 5), spikeTrain((spikeSpeed < 5))/Fs, '.','color', [0 191 255]/255,'markersize', 0.5)

hold on 
plot(spikePos(spike.lap > 0), spikeTrain((spike.lap > 0))/Fs, '.','color', 'r','markersize', 0.5)

ylabel('Time(Sec)', 'fontsize', 16)
xrange = xlim;

h1 = subplot(5,1,4);
[n1 pb] = hist(ripplePos, 50);
bar(pb, n1)
set(h1,'xlim', xrange)
ylabel('#Ripples', 'fontsize', 12)
yrange1 = n1;


h2 = subplot(5,1,5);
[n2 pb]= hist(eventPos, 50);
bar(pb, n2)
set(h2, 'xlim',xrange)
xlabel('Position(cm)', 'fontsize', 16)
ylabel('#PopBurst', 'fontsize', 12)
yrange2 = n2;

axis([h1 h2],[xrange 0 max([n1 n2])])

