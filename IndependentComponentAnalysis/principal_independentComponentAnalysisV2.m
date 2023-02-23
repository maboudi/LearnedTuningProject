% This script in intended to calcualte ICs and their corresponding weights
% from RUN period and their activation strengths during sleep periods: PRE and POST. 
% This analysis is based on two time scales: 250ms to detect reactivations of coarse
% time-scale assemblies and 20ms to detect fine time-scale firing regimes,
% each possiblty corresponding to a track position or in general a neuronal
% population state. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing spike-times; dividing each period to time bins with an
% amount of overlap between successive time bins

binDur = 0.250; % in sec
binOverlapRatio = 0;
stepSize = (1-binOverlapRatio)*binDur;

avgWinDur  = 15*60; % running average window in sec
avgWinSize = avgWinDur/binDur;
avgWindow  = ones(avgWinSize, 1)/avgWinSize;

posBinSize = 2; % in cm

binCenters   = struct('PRE', [], 'RUN', [], 'POST', []);
binnedFiring = struct('PRE', [], 'RUN', [], 'POST', []);

% Regarding the RUN Period, we are not exluding any part; active RUN and
% quiet awake both are included

% RUN
runSpeedThresh = 0;
binOverlapRatio = 0;
[runBinnedfiring, posbinIdx, linposcenters] = binRunData_1D(spikeStruct, behavior.time(2, :), behavior, speed, runSpeedThresh, binDur, binOverlapRatio, posBinSize, fileinfo);

% PRE
qclus = [1 2 3];
binCenters.PRE   = (behavior.time(1, 1):stepSize:(behavior.time(1, 2)-binDur)) + binDur/2;
binnedFiring.PRE = timeBinning_withOverlap(behavior.time(1, :) , spikeStruct, qclus, binDur, binOverlapRatio, fileinfo); % the theta periods

% RUN
binCenters.RUN   = (behavior.time(2, 1):stepSize:(behavior.time(2, 2)-binDur)) + binDur/2;
binnedFiring.RUN = timeBinning_withOverlap(behavior.time(2, :) , spikeStruct, qclus, binDur, binOverlapRatio, fileinfo); % the theta periods

% POST
binCenters.POST   = (behavior.time(3, 1):stepSize:(behavior.time(3, 2)-binDur)) + binDur/2;
binnedFiring.POST = timeBinning_withOverlap(behavior.time(3, :) , spikeStruct, qclus, binDur, binOverlapRatio, fileinfo); % the theta periods


% pooledPBEsRUN = [];
% for ii = 1: size(primaryPBEs_run, 1)
%     
%     firstlastbins = floor((primaryPBEs_run(ii, 1:2)-behavior.time(2,1))*1000/20);
%     pooledPBEsRUN = [pooledPBEsRUN; [firstlastbins(1):firstlastbins(2)]'];  
%     
% end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the RUN ICs and IC weights

[ICs, ICweights, eigenValues,thresh, mu, sigma] = pca_ica(runBinnedfiring); 
nICs = size(ICs, 1);

for ii = 1:nICs; ICweights(ii, :) = ICweights(ii, :)/norm(ICweights(ii, :)); end

ICs_sign = sum(ICweights, 2)./abs(sum(ICweights, 2));
ICweights = ICweights .* repmat(ICs_sign, [1, size(ICweights, 2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating the (re)activation strengths

activationStrength = struct('PRE', [], 'RUN', [], 'POST', []);
activationStrength = structfun(@(x) struct('original', [], 'smoothed', []), activationStrength, 'UniformOutPut', false);

% RUN
activationStrength.RUN.original  = calActStrength(ICweights, binnedFiring.RUN{1,2}(:, pooledPBEsRUN));
activationStrength.RUN.smoothed  = arrayConv(activationStrength.RUN.original, avgWindow); % running average

% PRE
activationStrength.PRE.original  = calActStrength(ICweights, binnedFiring.PRE{1,2}); 
activationStrength.PRE.smoothed  = arrayConv(activationStrength.PRE.original, avgWindow);

% POST
activationStrength.POST.original = calActStrength(ICweights, binnedFiring.POST{1,2}); 
activationStrength.POST.smoothed = arrayConv(activationStrength.POST.original, avgWindow);


% zscore the activation strength time course (the 15 minutes smoothed one)
% in respect to the PRE distribution

PREmean = mean(activationStrength.PRE.smoothed, 2);
PREstd  = std(activationStrength.PRE.smoothed, [], 2);

activationStrength.PRE.smoothed = zscore_KM(activationStrength.PRE.smoothed, PREmean, PREstd);
activationStrength.RUN.smoothed = zscore_KM(activationStrength.RUN.smoothed, PREmean, PREstd);
activationStrength.POST.smoothed = zscore_KM(activationStrength.POST.smoothed, PREmean, PREstd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sort the ICs based on their calculated spatial tunings

nPosBins = size(spatialTunings_LR, 2);
RUNbinsPostPr = struct('RL', [], 'LR', [], 'integrated', []);

spatialTunings_LR(spatialTunings_LR==0) = 1e-4;
RUNbinsPostPr.LR = baysDecoder(binnedFiring.RUN{1,2}(:, pooledPBEsRUN), spatialTunings_LR, binDur);

spatialTunings_RL(spatialTunings_RL==0) = 1e-4;
RUNbinsPostPr.RL = baysDecoder(binnedFiring.RUN{1,2}(:, pooledPBEsRUN), spatialTunings_RL, binDur);

RUNbinsPostPr.integrated = (RUNbinsPostPr.LR + RUNbinsPostPr.RL) ./ repmat(sum(RUNbinsPostPr.LR + RUNbinsPostPr.RL), [nPosBins, 1]);

RUNbinsPostPr.LR = RUNbinsPostPr.LR ./ repmat(sum(RUNbinsPostPr.LR, 1), [size(RUNbinsPostPr.LR, 1) 1]);
RUNbinsPostPr.RL = RUNbinsPostPr.RL ./ repmat(sum(RUNbinsPostPr.RL, 1), [size(RUNbinsPostPr.RL, 1) 1]);



ICAvgPosterior = zeros(nPosBins, nICs);
for ic = 1:nICs
   
    currWeights = activationStrength.RUN.original(ic, :);
    
    ICAvgPosterior(:, ic) = sum(repmat(currWeights, [nPosBins 1]) .* RUNbinsPostPr.integrated , 2) ./ sum(currWeights);

end

[~, ICpeakPositions] = max(ICAvgPosterior);
[~, ICpositionsortInd] = sort(ICpeakPositions, 'ascend');

figure; 
imagesc(ICAvgPosterior(:, ICpositionsortInd), [-0.4 0.4]); colormap('jet')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate spatial tunings of the IC (through reactivation strengths)
% what would be the best way to calculate the spatial tunings? (the ICs themselves or their activation strenghts)


uniqPositions = unique(posbinIdx(:, 1));
ICtunings     = zeros(nICs, numel(uniqPositions));

for ii = 1: numel(uniqPositions); ICtunings(:, ii) = sum(activationStrength.RUN.original(:, posbinIdx(:, 1) == uniqPositions(ii)), 2); end



%%%%%%%%%%
% figures

figure;
set(gcf, 'position', [2500 700 620 750])

% eigenValues
subplot(2,5,[1 2])

hold on
plot(1:numel(eigenValues), eigenValues, 'linewidth', 2, 'color', 'k')

signPCidx  = find(eigenValues > thresh);
insigPCidx = setdiff(1:numel(eigenValues), signPCidx); 

plot(signPCidx, eigenValues(signPCidx), 's', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', 3)
plot(insigPCidx, eigenValues(insigPCidx), 'o', 'markerfacecolor', 'none', 'markeredgecolor', 'k', 'markersize',1 )
line([0 numel(eigenValues)], [thresh thresh], 'linestyle','--', 'linewidth', 1, 'color','r')

hold off

set(gca, 'fontsize', 10, 'box', 'off', 'linewidth', 1)
ylabel('eigenvalues', 'fontsize', 12)
xlabel('principal components', 'fontsize', 12)


% ICs
subplot(2,5, 3:5)
cm = redblue;

imagesc(ICweights', [-max(abs(ICweights(:))) max(abs(ICweights(:)))]); colormap(cm); colorbar; 
set(gca, 'ydir', 'normal', 'fontsize', 10, 'linewidth', 1)
xlabel('independent components', 'fontsize', 12)
ylabel('units', 'fontsize', 12)


% sptial tuning of the ICs
subplot(2,5,8:10)

ICtunings_norm = zscore(ICtunings(:, 2:end), [], 2); % excluding the zero position bin and z scoring
imagesc(1:size(ICtunings_norm, 2)*2, 1:nICs,  ICtunings_norm); colormap(cm); colorbar;
set(gca, 'ydir', 'normal', 'fontsize', fontsize)
xlabel('track position(cm)', 'fontsize', 12)
ylabel('independent components', 'fontsize', 12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactivation strength of IC patterns

plotICactivationStrengths(activationStrength, binCenters, behavior, 8, 4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pooled reactivation strengths

% calculte the average reactivations corresponding to each conditions (PRE, RUN, and POST) by pooling over all the IC patterns

activationStrength2.PRE.original = mean(activationStrength.PRE.original, 1);
activationStrength2.RUN.original = mean(activationStrength.RUN.original, 1);
activationStrength2.POST.original = mean(activationStrength.POST.original, 1);

activationStrength2.PRE.smoothed = conv(activationStrength2.PRE.original, avgWindow, 'same');
activationStrength2.RUN.smoothed = conv(activationStrength2.RUN.original, avgWindow, 'same');
activationStrength2.POST.smoothed = conv(activationStrength2.POST.original, avgWindow, 'same');



% zscoring the pooled average reactivations in different conditions in respect to PRE distribution

PREmean = mean(activationStrength2.PRE.smoothed, 2);
PREstd  = std(activationStrength2.PRE.smoothed, [], 2);

activationStrength2.PRE.smoothed = zscore_KM(activationStrength2.PRE.smoothed, PREmean, PREstd);
activationStrength2.RUN.smoothed = zscore_KM(activationStrength2.RUN.smoothed, PREmean, PREstd);
activationStrength2.POST.smoothed = zscore_KM(activationStrength2.POST.smoothed, PREmean, PREstd); 

% plot the Pooled reactivation strength
plotICactivationStrengths(activationStrength2, binCenters, behavior, 1, 1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finding the lags between pairs of IC

nTimeBins = size(activationStrength.RUN.original, 2);
binIdx = (1:nTimeBins)';

jitterLength = 10;
nIter = 100; % number of shuffle distributions


shuffleBinIdx = repmat(binIdx, [1 nIter]);
tic

%%% calculating jitter data
for jj= 1:nIter
    jj
    alreadyUsedBins = [];
    sign = -1;
    
    for ii = 1: floor(nTimeBins/2)

        availableBins = setdiff(binIdx, alreadyUsedBins);

        bin1 = availableBins(randi(numel(availableBins)));

        flag = 0;


        if find(ismember(setdiff((bin1-jitterLength):(bin1+jitterLength), bin1), availableBins))
            while flag == 0
                sign = sign*-1;
                bin2 = bin1 + randi(jitterLength)*sign;

                if ismember(bin2, availableBins)
                   flag = 1;
                end
            end

            shuffleBinIdx(bin1, jj) = bin2;
            shuffleBinIdx(bin2, jj) = bin1;

            alreadyUsedBins = [alreadyUsedBins; bin1; bin2];
        else
            alreadyUsedBins = [alreadyUsedBins; bin1];
        end

    end
end

toc


nICs2plot = 10;
ICs2Plot = randperm(nICs, nICs2plot);

ICs2Plot = [1:10]; 

figure; 

for ii = 1:nICs2plot
    flag = 0;
    for jj = ii:nICs2plot
        

        % actual
        [xcorr_actual, lags] = xcorr(activationStrength.RUN.original(ICs2Plot(ii),binIdx), activationStrength.RUN.original(ICs2Plot(jj),binIdx));
        
        xcorr_null = zeros(nIter, length(xcorr_actual));
        for iIter = 1:nIter
           xcorr_null(iIter, :) = xcorr(activationStrength.RUN.original(ICs2Plot(ii), binIdx), activationStrength.RUN.original(ICs2Plot(jj), shuffleBinIdx(:, iIter))); 
        end
        
        
        subplot(nICs2plot, nICs2plot, (ii-1)*nICs2plot+jj)
        set(gca, 'fontsize', 5)
        hold on
        
        if jj == ii
            bar(lags, xcorr_actual, 'FaceColor', 'b', 'EdgeColor', 'b')
        else
            bar(lags, xcorr_actual, 'FaceColor', 'k', 'EdgeColor', 'k')
        end
        
        plot(lags, median(xcorr_null, 1), 'color', 'm', 'linewidth', 1)
        plot(lags, prctile(xcorr_null, 95, 1), 'color', 'm', 'linestyle', ':')
        plot(lags, prctile(xcorr_null, 5, 1), 'color', 'm', 'linestyle', ':')
        
        hold off
        
        xlim([-20 20])
       
        if ii == 1; title(sprintf('IC #%d', ICs2Plot(jj)), 'FontWeight', 'normal'); end
        
        if flag == 0 
           xlabel('time bin')
           ylabel({sprintf('IC #%d', ICs2Plot(ii));'';'cross-corr.'})
           flag =1;
        end
           
        xticks([-20 -15 -10 -5 0 5 10 15 20])
    end
end


%% finding the lags between IC pairs version 2
% This time using only the RUN PBEs to calculate cross-correlations



nPBEs_PRE = size(PREbinnedPBEs.data, 1);

% Data = POSTbinnedPBEs.data(randperm(nPBEs_POST, 1000), :);
% nPBEs_POST = 1000;
activationStrPBEs = cell(nPBEs_PRE, 1);
for pbe = 1:nPBEs_PRE
    activationStrPBEs{pbe} = calActStrength(ICweights, PREbinnedPBEs.data{pbe, 2});
end


maxLag    = 20; 
nShuffles = 10;

xcorr_actual_PBEsummed = cell(nICs2plot);
temporalBias           = zeros(nICs2plot);

xcorr_shuffle_PBEsummed = cell(nICs2plot);
temporalBias_shuffle    = cell(nICs2plot);


for ii = 1:nICs
    ii
    for jj = ii:nICs
  
        % actual
        xcorr_actual_PBE = zeros(nPBEs_PRE, 2*maxLag+1);

        for pbe = 1:nPBEs_PRE
            xcorr_actual_PBE(pbe, :) = xcorr(activationStrPBEs{pbe}(ii, :), activationStrPBEs{pbe}(jj, :), maxLag);
        end
        
        xcorr_actual_PBEsummed{ii, jj} = sum(xcorr_actual_PBE, 1);
        
        A = sum(xcorr_actual_PBEsummed{ii, jj}(1:maxLag));
        B = sum(xcorr_actual_PBEsummed{ii, jj}(maxLag+2:2*maxLag+1));
        temporalBias(ii, jj) = (A-B)/(A+B);
        
        
        % shuffle
        xcorr_shuffle_PBE = zeros(nShuffles, 2*maxLag+1, nPBEs_PRE);
        for pbe = 1:nPBEs_PRE
   
            s1 = activationStrPBEs{pbe}(ii, :);
            s2 = activationStrPBEs{pbe}(jj, :);
            
            nTimeBins = size(s1, 2);
            
            temp = zeros(nShuffles, 2*maxLag+1);
            
            for ishuffle = 1:nShuffles % why parfor is so slow with this??
                temp(ishuffle, :) = xcorr(s1, s2(randperm(nTimeBins)), maxLag);
            end
            
            xcorr_shuffle_PBE(:,:, pbe) = temp;    
        end
        
        xcorr_shuffle_PBEsummed{ii, jj} = sum(xcorr_shuffle_PBE, 3);
        
        temporalBias_shuffle{ii, jj}    = zeros(nShuffles, 1);
        
        for ishuffle = 1:nShuffles
            
            A = sum(xcorr_shuffle_PBEsummed{ii, jj}(ishuffle, 1:maxLag));
            B = sum(xcorr_shuffle_PBEsummed{ii, jj}(ishuffle, maxLag+2:2*maxLag+1));

            temporalBias_shuffle{ii, jj}(ishuffle) = (A-B)/(A+B);    
        end  

    end
end 


ICs2Plot  = 11:20;
% ICs2Plot  = randperm(nICs, 10);
nICs2Plot = numel(ICs2Plot);
lags = -maxLag:maxLag;


for ii = 1:nICs2Plot
    flag = 0;
   for jj = ii:nICs2Plot
      
        subplot(nICs2plot, nICs2plot, (ii-1)*nICs2plot+jj)
        set(gca, 'fontsize', 6)
        hold on
        
        if jj == ii
            bar(lags, xcorr_actual_PBEsummed{ICs2Plot(ii), ICs2Plot(jj)}, 'FaceColor', 'b', 'EdgeColor', 'b')
        else
            bar(lags, xcorr_actual_PBEsummed{ICs2Plot(ii), ICs2Plot(jj)}, 'FaceColor', 'k', 'EdgeColor', 'k')
        end
        
        plot(lags, median(xcorr_shuffle_PBEsummed{ICs2Plot(ii), ICs2Plot(jj)}, 1), 'color', 'm', 'linewidth', 1)
        plot(lags, prctile(xcorr_shuffle_PBEsummed{ICs2Plot(ii), ICs2Plot(jj)}, 95, 1), 'color', 'm', 'linestyle', ':')
        plot(lags, prctile(xcorr_shuffle_PBEsummed{ICs2Plot(ii), ICs2Plot(jj)}, 5, 1), 'color', 'm', 'linestyle', ':')

        
        hold off
        
        xlim([-20 20])
       
        if ii == 1; title(sprintf('IC #%d', ICs2Plot(jj)), 'FontWeight', 'normal'); end
        
        if flag == 0 
           xlabel('time bin')
           ylabel({sprintf('IC #%d', ICs2Plot(ii));'';'cross-corr.'})
           flag =1;
        end
           
        xticks([-20 -10  0  10  20])
        
        yLimit = ylim;
        
        pval = length(find(abs(temporalBias_shuffle{ICs2Plot(ii), ICs2Plot(jj)}) >= abs(temporalBias(ICs2Plot(ii), ICs2Plot(jj)))))/(nShuffles+1);
%         temporalBias_zscore = (temporalBias{ICs2Plot(ii), ICs2Plot(jj)} - mean(temporalBias_shuffle{ICs2Plot(ii), ICs2Plot(jj)}))/std(temporalBias_shuffle{ICs2Plot(ii), ICs2Plot(jj)});
        text(5, yLimit(1)+0.8*range(yLimit), {sprintf('tb= %.2f', temporalBias(ICs2Plot(ii), ICs2Plot(jj))); sprintf('pval= %0.2f', pval)}, 'fontsize', 7)
       
   end
    
end

save('PREtemporalBias_withoutNormalization.mat', 'xcorr_actual_PBEsummed', 'xcorr_shuffle_PBEsummed', 'temporalBias', 'temporalBias_shuffle')


%% plotting figures related to temporal biases for each period and correlation between them

figure; 
set(gcf, 'position', [100 100 900 260])

colorLimit = [-1 1];

subplot(2,3,1)
imagesc(PREtemporalBias.temporalBias, [-1 1]); colormap(cm)
set(gca, 'linewidth', 2)
title({'PRE'; '';'IC j'}, 'fontsize', 12, 'fontweight', 'normal')
ylabel('IC i', 'fontsize', 12)

h2 = subplot(2,3,2);
imagesc(RUNtemporalBias.temporalBias, [-1 1]); colormap(cm)
set(gca, 'linewidth', 2)
title({'RUN'; '';'IC j'}, 'fontsize', 12, 'fontweight', 'normal')
ylabel('IC i', 'fontsize', 12)



subplot(2,3,3);
imagesc(POSTtemporalBias.temporalBias, [-1 1]); colormap(cm)

set(gca, 'linewidth', 2)
title({'POST'; '';'IC j'}, 'fontsize', 12, 'fontweight', 'normal')
ylabel('IC i', 'fontsize', 12)

h = colorbar;
h.Position(1) = h.Position(1) + 0.08; 
h


RUN_PRE_tbCorr  = zeros(nICs, 1);
RUN_POST_tbCorr = zeros(nICs, 1);

RUN_POST_tbCorr_shuffle = zeros(nICs, 1000);
RUN_PRE_tbCorr_shuffle  = zeros(nICs, 1000);

for ii = 1: nICs
    
    PREtbs = [PREtemporalBias.temporalBias(ii, min(nICs, ii+1): nICs) PREtemporalBias.temporalBias(1:ii-1, ii)'];
%     PREtbs = PREtbs(PREtbs ~= 0);
    PREtbs = PREtbs./abs(PREtbs);
    
    RUNtbs = [RUNtemporalBias.temporalBias(ii, min(nICs, ii+1): nICs) RUNtemporalBias.temporalBias(1:ii-1, ii)'];
%     RUNtbs = RUNtbs(RUNtbs ~= 0);
    RUNtbs = RUNtbs./abs(RUNtbs);
    
    POSTtbs = [POSTtemporalBias.temporalBias(ii, min(nICs, ii+1): nICs) POSTtemporalBias.temporalBias(1:ii-1, ii)'];
%     POSTtbs = POSTtbs(POSTtbs ~= 0);
    POSTtbs = POSTtbs./abs(POSTtbs);
    
    
    RUN_PRE_tbCorr(ii)  = (PREtbs  * RUNtbs')/norm(PREtbs) /norm(RUNtbs);
    RUN_POST_tbCorr(ii) = (POSTtbs * RUNtbs')/norm(POSTtbs)/norm(RUNtbs);
    
    for ishuffle = 1:1000
        
        rndgen = randperm(numel(RUNtbs));
        RUN_POST_tbCorr_shuffle(ii, ishuffle) = (POSTtbs * RUNtbs(rndgen)')/norm(POSTtbs)/norm(RUNtbs(rndgen));  
        RUN_PRE_tbCorr_shuffle(ii, ishuffle)  = (PREtbs * RUNtbs(rndgen)')/norm(PREtbs)/norm(RUNtbs(rndgen));
    end
    
    
        
end


subplot(2,3,4)

pooled = [RUN_PRE_tbCorr; RUN_POST_tbCorr; RUN_PRE_tbCorr_shuffle(:); RUN_POST_tbCorr_shuffle(:)];
bins = linspace(min(pooled), max(pooled), 20);

hPRE = hist(RUN_PRE_tbCorr, bins);
hPRE = hPRE./sum(hPRE);
hPOST        = hist(RUN_POST_tbCorr, bins);
hPOST= hPOST./sum(hPOST); 

hPRE_shuffle = hist(RUN_PRE_tbCorr_shuffle(:), bins);
hPRE_shuffle = hPRE_shuffle./sum(hPRE_shuffle);
hPOST_shuffle = hist(RUN_POST_tbCorr_shuffle(:), bins);
hPOST_shuffle = hPOST_shuffle./sum(hPOST_shuffle);


bar(bins, hPRE, 'FaceColor', [255 233 64]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
hold on
bar(bins, hPOST, 'FaceColor', [64 86 255]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
xlabel({'correlation'; 'with RUN temporal biases'}, 'fontsize', 12)
ylabel('probability', 'fontsize', 12)
legend('RUN-PRE', 'RUN-POST', 'location', 'northwest', 'box', 'off')
set(gca, 'linewidth', 2, 'box', 'off')



subplot(2,3,5)

bar(bins, hPRE_shuffle, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.6)
hold on
bar(bins, hPRE, 'FaceColor', [255 233 64]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
xlabel({'correlation'; 'with RUN temporal biases'}, 'fontsize', 12)
ylabel('probability', 'fontsize', 12)
title('RUN-PRE', 'fontsize', 12)
legend('shuffle', 'data', 'location', 'northwest', 'box', 'off')
set(gca, 'linewidth', 2, 'box', 'off')



subplot(2,3,6)

bar(bins, hPOST_shuffle, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.6)
hold on
bar(bins, hPOST, 'FaceColor', [64 86 255]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
xlabel({'correlation'; 'with RUN temporal biases'}, 'fontsize', 12)
ylabel('probability', 'fontsize', 12)
title('RUN-POST', 'fontsize', 12)
legend('shuffle', 'data', 'location', 'northwest', 'box', 'off')
set(gca, 'linewidth', 2, 'box', 'off')


PREtemporalBiasAll = [];
RUNtemporalBiasAll = [];
POSTtemporalBiasAll = [];

for ii = 1:nICs-1
    PREtemporalBiasAll = [PREtemporalBiasAll; diag(PREtemporalBias.temporalBias, ii)];
    RUNtemporalBiasAll = [RUNtemporalBiasAll; diag(RUNtemporalBias.temporalBias, ii)];
    POSTtemporalBiasAll = [POSTtemporalBiasAll; diag(POSTtemporalBias.temporalBias, ii)];
end


figure; 

subplot(1,2,1)
set(gca, 'linewidth', 2, 'box', 'off')
plot(RUNtemporalBiasAll, PREtemporalBiasAll, '.k', 'markersize', 3)
x = [ones(length(RUNtemporalBiasAll), 1) RUNtemporalBiasAll];
b= x\PREtemporalBiasAll;
ycalc = x*b;
hold on
plot(RUNtemporalBiasAll, ycalc, 'color', [255 233 64]/255, 'linewidth',2)
xlabel('RUN temporal bias', 'fontsize', 12)
ylabel('PRE temporal bias', 'fontsize', 12)
xlim([-1 1])
ylim([-20 20])

subplot(1,2,2)
set(gca, 'linewidth', 2, 'box', 'off')
plot(RUNtemporalBiasAll, POSTtemporalBiasAll, '.', 'markersize', 3)
b= x\POSTtemporalBiasAll;
ycalc = x*b;
hold on
plot(RUNtemporalBiasAll, ycalc, 'color', [64 86 255]/255, 'linewidth',2)
xlabel('RUN temporal bias', 'fontsize', 12)
ylabel('POST temporal bias', 'fontsize', 12)
xlim([-1 1])
ylim([-1 1])


%% visualizing how population firings, posterior probability and reactivation strength of different ICs correlate with each other

% activation strength during the PBEs
nPBEs = length(POSTbinnedPBEs.data(:, 2));
activationStrPBEs = cell(nPBEs, 1);

for ii = 1:nPBEs
   
    activationStrPBEs{ii} = calActStrength(ICweights, POSTbinnedPBEs.data{ii, 2});
end


BDscores = BDseqscore.POST.data.unitIDshuffle.replayScore.prctilescore;

% PBE2Plot = randperm(numel(BDscores), 40);

PBE2Plot = find(BDscores > 95);
PBE2Plot = PBE2Plot(randperm(numel(PBE2Plot), 40));

%%

period = 'POST';
filename = [period '_ICsActivationDuringPBEs_highBDscores2'];
directory = fullfile(mainDir, 'ICA', period);
mkdir(directory)

textOnTop = sprintf('Session: %s\nPeriod: %s', sessionName, period);

plotICactivation4PBEs(POSTbinnedPBEs.data, activationStrPBEs, posteriorProbMatrix.POST.data, BDscores, runTemplate_LR, runTemplate_RL, PBE2Plot, directory, filename, textOnTop)


%%


%%%%%%%%%%%
% functions


function activationStrength = calActStrength(ICweights, x)

x = zscore(x, [], 2);
nTimebins = size(x,2);
nICs = size(ICweights, 1); % rows correspond to weights for each IC

activationStrength = zeros(nICs, nTimebins);

for ii = 1:nICs
    currentICbase = ICweights(ii, :)';
    ProjectorOp = currentICbase * currentICbase';
    ProjectorOp = ProjectorOp - ProjectorOp.*eye(size(ProjectorOp)); % removing the diagonal to prevent high activation strength casued by the isolated activity of a single neuron with high weight in the pattern
    
    for jj = 1:nTimebins
        activationStrength(ii, jj) = x(:, jj)' * ProjectorOp * x(:, jj);
    end
end

end

function smoothedMat = arrayConv(inputMat, avgWindow)

smoothedMat = zeros(size(inputMat));
for ii = 1: size(inputMat, 1)
    smoothedMat(ii, :) = conv(inputMat(ii, :), avgWindow, 'same');
end

end


function z = zscore_KM(x, mu, sigma)

z = (x - repmat(mu, [1 size(x, 2)]))./repmat(sigma, [1 size(x, 2)]);

end
