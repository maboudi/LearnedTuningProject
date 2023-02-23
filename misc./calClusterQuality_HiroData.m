function  clusterQuality = calClusterQuality(baseFile)

% mainDir = '/data/GrossmarkData_BapunsClusters';
% % mainDir = '/data/Hiro_clusters/Roy_Data/Roy-maze3';
% currentSession = 'Buddy_06272013';
% 
% datasetDir = fullfile(mainDir, currentSession);
% cd(datasetDir)
% allFet = dir([currentSession '.fet.*']);

[filePath, name] = fileparts(baseFile);

fileNameBase = fullfile(filePath, name);

try
    Par = LoadXml([fileNameBase '.xml']);
catch
    Par = LoadXml(fileNameBase);
end


if isfield(Par, 'SpkGrps')
    shankInfo = Par.SpkGrps;
elseif isfield(Par, 'ElecGp')
    shankInfo = Par.ElecGp;
elseif isfield(Par, 'AnatGrps')
    shankInfo = Par.AnatGrps;
else
    error('\nShank info was not found!')
end    
nShanks = numel(shankInfo);



clusterQuality = struct('mahalDist', [], 'isoDist', [], 'Lratio', [], 'clus', []);



for iShank = 1:nShanks
    
    iShank
    if isfield(shankInfo, 'Skip')
        channs = shankInfo(iShank).Channels(~shankInfo(iShank).Skip);
    else
        channs = shankInfo(iShank).Channels;
    end
%     channs = shankInfo(iShank).Channels;
%     channs(ismember(channs, badChannels)) = [];
    
    nChannels = numel(channs);    
    
    
    try  
        fileID = fopen([fileNameBase '.fet.' num2str(iShank)]);
        featureFile = fscanf(fileID, '%d');
    catch
        fprintf('features file is missing for shank %d!', iShank)
        continue
    end
    
    
    nFeatures    = featureFile(1);
    featureFile2 = featureFile(2:end);
    
    try
        featureFile2 = reshape(featureFile2, [nFeatures numel(featureFile2)/nFeatures]);
    catch
        fprintf('\nCheck shank %d something wrong with the number of channels or features!', iShank) 
        continue
    end
    

    try
        featureFile2 = featureFile2(1:nChannels*3, :);
    catch
        fprintf('\nCheck shank %d', iShank)
        featureFile2 = featureFile2(1:nFeatures, :);
    end
    
    
%     resFile = load([baseFile '.res.' num2str(iShank)]);
    
%     spkFileName = [fileNameBase '.spk.' num2str(iShank)];
%     
%     if exist(spkFileName, 'file')
%         [spikeWaveForms, nChannels] = readSpk(spkFileName, nChannels, 48); % spikeWaveForms = [nChannels * nSamplesPerSpike * nSpikes]
%     else
%         continue
%     end
%     
%     ifanyAmp = squeeze(sum(sum(spikeWaveForms, 1), 2));
%     spikeWaveForms = spikeWaveForms(:, :, ifanyAmp~=0);
%     
%     
%     nSpikes = size(spikeWaveForms, 3);
%     
%     % calculating the PCs, power and the range of the 
%     
%     PCAs_global = zeros(3, nChannels, nSpikes);
%     for k = 1:nChannels % for each channel
%          [~, scores] = pca(permute(spikeWaveForms(k,:,:),[3,2,1]),'NumComponents',3); % do we need the z-scoring here before PCA??
%          PCAs_global(:,k,:) = scores';
%     end
%     featureFile2 = reshape(PCAs_global, [nChannels*3, nSpikes]);
%     
    
    cluFile = load([fileNameBase '.clu.' num2str(iShank)]);
    
    nClus = cluFile(1);
    cluFile = cluFile(2:end);
    cluFile = cluFile(ifanyAmp ~= 0);
    uniqClus = unique(cluFile);
    

    cluFeatures   = cell(nClus, 1);
% %     spkWv         = cell(nClus, 1);
%     spkWv_mean    = cell(nClus, 1);
%     spkWv_peakAmp = cell(nClus, 1);
    
%     nClusterSpikes = zeros(nClus, 1);
    for iClu = 1:nClus
        currClu = uniqClus(iClu);
        cluFeatures{iClu}    = featureFile2(:, cluFile == currClu);
%         nClusterSpikes(iClu) = numel(find(cluFile == currClu));
        
        
% %         currSpkWv = spikeWaveForms(:,:, cluFile == currClu);
%         
% %         spkWv{iClu} = currSpkWv;
%         spkWv_mean{iClu} = nanmean(spikeWaveForms(:,:, cluFile == currClu), 3);
%         spkWv_peakAmp{iClu} = min(spkWv_mean{iClu}, [], 2);
        
    end
    
%     clear spikeWaveForms
% %     clear spkWv 
    
    df = nChannels*3;
    clusterQuality(iShank).mahalDist = cell(nClus, nClus);
    clusterQuality(iShank).Lratio    = nan(nClus, nClus);
    clusterQuality(iShank).isoDist   = nan(nClus, nClus);
    
    
    for iClu = 1:nClus
        
       for jClu = setdiff(1:nClus, iClu)
           
           
           % find up to 4 channels with maximum difference between the two
           % units
           
           nSelectChann = 4;
           
           spkAmpDiffs = abs(spkWv_peakAmp{jClu} - spkWv_peakAmp{iClu});
           [~, sortIdx] = sort(spkAmpDiffs, 'descend');
           channsWithMaxDiffs = sortIdx(1:nSelectChann);
           
           
           featuresIdx2Include = nan(nSelectChann*3, 1);
           for ch = 1:nSelectChann
               firstFeature = (channsWithMaxDiffs(ch)-1)*3+1;
               featuresIdx2Include((ch-1)*3+1:ch*3) = firstFeature:(firstFeature+2);
           end
           df = nSelectChann*3;
           
           refClu = cluFeatures{iClu}(featuresIdx2Include, :)';
           testClu = cluFeatures{jClu}(featuresIdx2Include, :)';
           
           
           try % in case there are not enough spikes (less than the number of features, it gives an error)
               
               allMahalDist = mahal(testClu, refClu);               
                
               L = sum(1-chi2cdf(allMahalDist, df));
               clusterQuality(iShank).Lratio(jClu, iClu) = L/nClusterSpikes(iClu);
               
%                L/nClusterSpikes(iClu)
               
  
               if nClusterSpikes(iClu) < nClusterSpikes(jClu)  
                    sorted = sort(clusterQuality(iShank).mahalDist{jClu, iClu});
                    clusterQuality(iShank).isoDist(jClu, iClu) = sorted(nClusterSpikes(iClu));
               elseif nClusterSpikes(iClu) > nClusterSpikes(jClu)
                    clusterQuality(iShank).isoDist(jClu, iClu) = median(clusterQuality(iShank).mahalDist{jClu, iClu});
               end
           catch
               clusterQuality(iShank).mahalDist{jClu, iClu} = nan;
               clusterQuality(iShank).Lratio(jClu, iClu)    = nan;
               clusterQuality(iShank).isoDist(jClu, iClu)   = nan;
           end
           
       end
        
    end
    clusterQuality(iShank).clus = uniqClus;  

%     clear PCAs_global cluFeatures cluFile featureFile2 ifanyAmp scores 
    
    
    
    
%     %% temp figure
%     
% %     iClu = 26;
% %     noiseClus = [2 31 27 22 10 20];
%     
%     iClu = 14;
%     noiseClus = setdiff(1:20, [iClu 12 16 5]);
%     
% 
% 
%     colors = colormap('lines');
%     
%     colors = colors(floor(linspace(1,64, numel(noiseClus))), :);
%     
% %     colors = colors([1 2 3 4 5 6 7 8 9 10], :);
%     
%   
% 
%     plotwidth  = 400;
%     plotheight = 350;
%     fontsize   = 6;
% 
%     f= figure;
%     clf(f);
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [plotwidth plotheight]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
%     
%     
%     nor = 3;
%     noc = 6;
% 
%     leftmargin = 40;  rightmargin = 40;    topmargin = 40;    bottommargin = 40;
% 
%     gapc = 15;
%     gapr = 10;
% 
%     sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);
%     
%     
%     tt = 0;
%     for jClu = noiseClus
%         
%         tt = tt+1;
%             
%         axes('position', sub_pos{tt},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
% 
%         hold on
%         
%         
%         spkAmpDiffs = abs(spkWv_peakAmp{jClu} - spkWv_peakAmp{iClu});
%         [~, sortIdx] = sort(spkAmpDiffs, 'descend');
%         channsWithMaxDiffs = sortIdx(1:nSelectChann);
% 
%         currChs = channsWithMaxDiffs(1:2);
% 
%         
%         
%         h1 = scatter(cluFeatures{iClu}((currChs(1)-1)*3+1, 1:3000), cluFeatures{iClu}((currChs(2)-1)*3+1, 1:3000), 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('unit %d', iClu));
%         h2 = scatter(cluFeatures{jClu}((currChs(1)-1)*3+1, 1:3000), cluFeatures{jClu}((currChs(2)-1)*3+1, 1:3000), 0.5, 'filled', 'MarkerFaceColor', colors(tt, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3, 'DisplayName', sprintf('unit %d', jClu));
% 
%         xlabel(sprintf('ch %d', currChs(1)), 'fontsize', fontsize)
%         ylabel(sprintf('ch %d', currChs(2)), 'fontsize', fontsize)
%         
%         
%         ch1Cloud = [cluFeatures{iClu}((currChs(1)-1)*3+1, 1:3000) cluFeatures{jClu}((currChs(1)-1)*3+1, 1:3000)];
%         ch2Cloud = [cluFeatures{iClu}((currChs(2)-1)*3+1, 1:3000) cluFeatures{jClu}((currChs(2)-1)*3+1, 1:3000)];
%         
%         xlim([prctile(ch1Cloud, 1) prctile(ch1Cloud, 99)])
%         ylim([prctile(ch2Cloud, 1) prctile(ch2Cloud, 99)])
%         
%         xl = xlim;
%         yl = ylim;
%         
% 
%         currLratio = clusterQuality(iShank).Lratio(jClu, iClu);
%         
%         if currLratio == 0
%             txt = 'L-ratio = 0';
%         elseif currLratio < 1e-3 
%             txt = sprintf('L-ratio = %.1e', currLratio);
%         else
%             txt = sprintf('L-ratio = %.3f', currLratio);
%         end
% 
% 
%         text(xl(1), yl(2)+0.25*range(yl), txt, 'fontsize', fontsize)
% 
%         set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'XTick', 0, 'YTick', 0)
%         
%         if tt == 1
%             h = legend([h1, h2], 'location', 'northoutside', 'box', 'off');
%             h.Position(2) = h.Position(2) + 0.045;
%         else
%             h = legend(h2, 'location', 'northoutside', 'box', 'off');
%             h.Position(2) = h.Position(2) + 0.02;
%         end
% 
%         
% 
%         axis square
%         grid on
% 
%     end
% 
% 
%     
%     
%     %% for RatU
%     
%     
%     electrodePositions = [...
%     -7.5    105    0; ...
%     7.5     97.5   0; ...
%     -7.5    90     0; ...
%     7.5     82.5   0; ...
%     -7.5    75     0; ...
%     7.5     67.5   0; ...
%     -7.5    60     0; ...
%     7.5     52.5   0; ...
%     -7.5    45     0; ...
%     7.5     37.5   0; ...
%     -7.5    30     0; ...
%     7.5     22.5   0; ...
%     -7.5    15     0; ...
%     7.5     7.5    0; ...
%     -7.5    0      0; ...
%     7.5     -7.5   0; ...
%     ];
% 
%     p_base_shank = [9 -22; 14 -10; 14 130; -14 130; -14 -2.5; 9 -22];
%     p_base_elec  = [-3.6 -3.6; 3.6 -3.6; 3.6 3.6; -3.6 3.6; -3.6 -3.6];
%     
% 
%     
%     plotwidth  = 400;
%     plotheight = 350;
%     fontsize   = 6;
% 
%     colors = colormap('lines');
%     
%     colors = colors(floor(linspace(1,64, numel(noiseClus))), :);
%     
% %     colors = colors([1 2 3 4 5 6 7 8 9 10], :);
%     
%     f= figure;
%     clf(f);
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [plotwidth plotheight]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
%     
%     nor = 3;
%     noc = 6;
% 
%     leftmargin = 40;  rightmargin = 40;    topmargin = 30;    bottommargin = 30;
% 
%     gapc = 10;
%     gapr = 0;
% 
%     sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);
% 
%     tt=0;
%     for jClu = noiseClus
%         
%         tt = tt + 1;
%         axes('position', sub_pos{tt},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
% 
%                 
%         hold on
%         p = p_base_shank;
%         pgon = polyshape(p);
% 
%         plot(pgon, 'FaceColor', [0.3 0.3 .3], 'EdgeColor', 'none')
% 
% 
%         for ielec = 1:nChannels
% 
%             p = p_base_elec + repmat(electrodePositions(ielec, 1:2), [5, 1]);
%             pgon = polyshape(p);
% 
%             plot(pgon, 'FaceColor', [1 1 1], 'EdgeColor', 'k')
%         end
%         
%         spkWv_iClu = spkWv_mean{iClu};
%         spkWv_jClu = spkWv_mean{jClu};
%         
%         maxApm = -min(spkWv_iClu(:)); % max amplitude across all channels on the shank
%         
%         for jj = 1:nChannels
% 
%             elecP = electrodePositions(jj, :);
% 
%            if mod(jj, 2) == 1
%                plot(elecP(1)-50-23: elecP(1)-50+24, spkWv_iClu(jj, :)*10/maxApm+elecP(2), 'color', [0 0 0 0.7], 'linewidth', 0.75)
%                plot(elecP(1)-2*50-23: elecP(1)-2*50+24, spkWv_jClu(jj, :)*10/maxApm+elecP(2), 'color', [colors(tt, :) 0.7], 'linewidth', 0.75)
% 
%            else
%                plot(elecP(1)+50-23: elecP(1)+50+24, spkWv_iClu(jj, :)*10/maxApm+elecP(2), 'color', [0 0 0 0.7], 'linewidth', 0.75)
%                plot(elecP(1)+2*50-23: elecP(1)+2*50+24, spkWv_jClu(jj, :)*10/maxApm+elecP(2), 'color', [colors(tt, :) 0.7], 'linewidth', 0.75)
% 
%            end
% 
%         end
%         axis tight
%          axis off
% 
%     end

    
end



% save(fullfile(datasetDir, [baseFile '.mahalDist.mat']), )
clusterQuality = rmfield(clusterQuality, 'mahalDist');
save(fullfile(filePath, [name '.clusterQuality.mat']), 'clusterQuality')



