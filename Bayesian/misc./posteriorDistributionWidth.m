
runPBEIdx = find(ismember({PBEInfo_replayScores.epoch}, {'run'}));
nPBEs = numel(runPBEIdx);

posteriorProbMat = cell(1, nPBEs);

for ipbe = 1:nPBEs
   
    posteriorProbMat{ipbe} = PBEInfo_replayScores(runPBEIdx(ipbe)).posteriorProbMat;
    
end

posteriorProbMat = cell2mat(posteriorProbMat);


posteriorProbMat(:, ~sum(posteriorProbMat, 1)) = [];


smoothWindow = gausswindow(3, 9);

nTimeBins = size(posteriorProbMat, 2);

peak   = zeros(nTimeBins, 1);
offset = zeros(nTimeBins, 1);
onset  = zeros(nTimeBins, 1);
posterior_width = zeros(nTimeBins, 1);
smoothed_PP = cell(nTimeBins, 1);

for ibin = 1:nTimeBins
    
    currPP = posteriorProbMat(:, ibin);
    currPP = conv(currPP, smoothWindow, 'same');
    
    [peaks, peakPos] = findpeaks(currPP);
    
    [peak(ibin), maxIdx] = max(peaks);
    peakPos = peakPos(maxIdx);
    
    
    % onset
    tempIdx = find(currPP(1:peakPos-1) <= 0.1*peak(ibin), 1, 'last');
    
    if ~isempty(tempIdx)
        onset(ibin) = tempIdx;
    else
        onset(ibin) = 1;
    end
    
    
    % offset
    tempIdx = find(currPP(peakPos+1: end) <= 0.1*peak(ibin), 1, 'first');
    if ~isempty(tempIdx)
        offset(ibin) = peakPos + tempIdx;
    else
        offset(ibin) = length(currPP);
    end
    
    
    posterior_width(ibin) = offset(ibin) - onset(ibin) + 1; 
    
    smoothed_PP{ibin} = currPP;
    
end

posterior_width(posterior_width > 0.9* size(posteriorProbMat, 1)) = nan;


figure; 
for ibin = 1:64
    
    subplot(8, 8, ibin)
    hold on

    plot(smoothed_PP{ibin}, 'linewidth', 2)
    yl = ylim;

    line([onset(ibin) onset(ibin)], yl, 'linewidth', 1, 'color', 'r')
    line([offset(ibin) offset(ibin)], yl, 'linewidth', 1, 'color', 'r')
    
    xticks([])
    yticks([])
    
    if ibin == 57
        xlabel('position')
        ylabel('posterior')
    end
    
end

suptitle('example maze PBE time bins')


figure; 
[a, b] = hist(posterior_width, 20);
bar(b*2,a, 'facecolor', 'k', 'edgecolor', 'none')

yl = ylim;
medianWidth = nanmedian(posterior_width)*2;

line([medianWidth medianWidth], yl, 'color', 'r', 'linewidth', 3)

xlabel('posterior width (cm)', 'fontsize', 14)
ylabel('# time bins', 'fontsize', 14)


%% 

figure; 


pbe = 292;

currPPM = PBEInfo_replayScores(pbe).posteriorProbMat;
begPos  = PBEInfo_replayScores(pbe).begPosition;
endPos  = PBEInfo_replayScores(pbe).endPosition; 


cl_min = min(currPPM(:));
cl_max = 0.5* max(currPPM(:));



subplot(1,2,1)
imagesc(currPPM, [cl_min cl_max])
colormap('hot')

line([begPos(1) endPos(1)], [begPos(2) endPos(2)], 'color', [1 1 1], 'linewidth', 2)
line([begPos(1) endPos(1)], [begPos(2) endPos(2)]-7, 'color', [1 1 1], 'linewidth', 1, 'lineStyle', ':')
line([begPos(1) endPos(1)], [begPos(2) endPos(2)]+7, 'color', [1 1 1], 'linewidth', 1, 'lineStyle', ':')

title('30 cm vicinity band')



subplot(1,2,2)
imagesc(currPPM, [cl_min cl_max])
colormap('hot')

line([begPos(1) endPos(1)], [begPos(2) endPos(2)], 'color', [1 1 1], 'linewidth', 2)
line([begPos(1) endPos(1)], [begPos(2) endPos(2)]-12, 'color', [1 1 1], 'linewidth', 1, 'lineStyle', ':')
line([begPos(1) endPos(1)], [begPos(2) endPos(2)]+12, 'color', [1 1 1], 'linewidth', 1, 'lineStyle', ':')

title('50 cm vicinity band')