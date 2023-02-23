clear;
clc;


periodNames = {'PRE'; 'RUN'; 'POST'};
shuffleMethods   = {'unitIDshuffle'; 'wPBEtimeswap'};
seqMetricMethods = {'weightedCorr'; 'replayScore'};



popBDscores = struct('PRE', [], 'RUN', [], 'POST', []);
for iperiod = 1:3
    
    currPeriod = (periodNames{iperiod});
    popBDscores.(currPeriod) = struct('unitIDshuffle', [], 'wPBEtimeswap', []);
    
    popBDscores.(currPeriod).unitIDshuffle = struct('weightedCorr', [], 'replayScore', []);
    popBDscores.(currPeriod).wPBEtimeswap = struct('weightedCorr', [], 'replayScore', []);
    
end


%% Grosmark


cd('/home/kouroshmaboudi/Documents/HMM_project/GreatLakes_firstRun_Nov2020/Grosmark_originalClusters')
currDir = pwd;

rr = dir;

for sessionNumber = 1:8
    sessionName = rr(sessionNumber+2).name;
    
    currSessionScores = fullfile(currDir, sessionName, 'BayesianDecoding');
    cd(currSessionScores)
    
    ll = dir;
    
    load(ll(3).name)
    

    datasets.PRE = PREbinnedPBEs;
    datasets.RUN = RUNbinnedPBEs;
    datasets.POST = POSTbinnedPBEs;

    periods = {'PRE';'RUN'; 'POST'};

    PBELen = struct('PRE', [], 'RUN', [], 'POST', []);

    for pp = 1:3

        currentPBEs = datasets.(periods{pp});
        nPBEs = size(currentPBEs, 1);

        PBELen.(periods{pp}) = zeros(nPBEs, 1);
        nFiringUnits.(periods{pp}) = zeros(nPBEs, 1);

        for pbe = 1:nPBEs

            PBELen.(periods{pp})(pbe) = size(currentPBEs{pbe, 2}, 2);

            nFiringUnits.(periods{pp})(pbe) = numel(find(sum(currentPBEs{pbe, 2}, 2)));

        end

    end


    secondaryPBEs_1.PRE = secondaryPBEs(secondaryPBEs(:, 1) > behavior.PREEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.PREEpoch(1, 2), :);
    secondaryPBEs_1.RUN = secondaryPBEs(secondaryPBEs(:, 1) > behavior.MazeEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.MazeEpoch(1, 2), :);
    secondaryPBEs_1.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.POSTEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.POSTEpoch(1, 2), :);


    acceptedEvts.PRE = find(secondaryPBEs_1.PRE(:, 6) == 1 & PBELen.PRE >= 4 & nFiringUnits.PRE >= 5);
    acceptedEvts.RUN = find(secondaryPBEs_1.RUN(:, 5) == 1 & PBELen.RUN >= 4 & nFiringUnits.RUN >= 5);
    acceptedEvts.POST = find(secondaryPBEs_1.POST(:, 6) == 1 & PBELen.POST >= 4 & nFiringUnits.POST >= 5);

    
    popBDscores_bySession.(sessionName) = BDseqscore;
    
    
    for ishuffleMethod = 1:2
        
        currShuffleMethod = shuffleMethods{ishuffleMethod};
        
        for iseqMetric = 1:2
            
            currSeqMetric = seqMetricMethods{iseqMetric};
            
            for iperiod = 1:3
                popBDscores.(periodNames{iperiod}).(currShuffleMethod).(currSeqMetric) = [popBDscores.(periodNames{iperiod}).(currShuffleMethod).(currSeqMetric); BDseqscore.(periodNames{iperiod}).data.(currShuffleMethod).(currSeqMetric).prctilescore(acceptedEvts.(periodNames{iperiod}))];  
            end
        end
    end
end


%% Hiro


cd('/home/kouroshmaboudi/Documents/HMM_project/GreatLakes_firstRun_Nov2020/Hiro')
currDir = pwd;

rr = dir;


for sessionNumber = 1:6
    
    
    sessionName = rr(sessionNumber+2).name;
    
    currSessionScores = fullfile(currDir, sessionName, 'BayesianDecoding');
    cd(currSessionScores)
    
    ll = dir;
    
    load(ll(3).name)
    
    
    datasets.PRE = PREbinnedPBEs;
    datasets.RUN = RUNbinnedPBEs;
    datasets.POST = POSTbinnedPBEs;

    periods = {'PRE';'RUN'; 'POST'};

    PBELen = struct('PRE', [], 'RUN', [], 'POST', []);

    for pp = 1:3

        currentPBEs = datasets.(periods{pp});
        nPBEs = size(currentPBEs, 1);

        PBELen.(periods{pp}) = zeros(nPBEs, 1);
        nFiringUnits.(periods{pp}) = zeros(nPBEs, 1);

        for pbe = 1:nPBEs

            PBELen.(periods{pp})(pbe) = size(currentPBEs{pbe, 2}, 2);

            nFiringUnits.(periods{pp})(pbe) = numel(find(sum(currentPBEs{pbe, 2}, 2)));

        end

    end
    

    secondaryPBEs_1.PRE = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(1, 1) & secondaryPBEs(:, 2) < behavior.time(1, 2), :);
    secondaryPBEs_1.RUN = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(2, 1) & secondaryPBEs(:, 2) < behavior.time(2, 2), :);
    secondaryPBEs_1.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(3, 1) & secondaryPBEs(:, 2) < behavior.time(3, 2), :);


    acceptedEvts.PRE = find(secondaryPBEs_1.PRE(:, 6) == 1 & PBELen.PRE >= 4 & nFiringUnits.PRE >= 5);
    acceptedEvts.RUN = find(secondaryPBEs_1.RUN(:, 5) == 1 & PBELen.RUN >= 4 & nFiringUnits.RUN >= 5);
    acceptedEvts.POST = find(secondaryPBEs_1.POST(:, 6) == 1 & PBELen.POST >= 4 & nFiringUnits.POST >= 5);

    
    sessionName(sessionName == '-') = '_';
    popBDscores_bySession.(sessionName) = BDseqscore;
    
    
    for ishuffleMethod = 1:2
        
        currShuffleMethod = shuffleMethods{ishuffleMethod};
        
        for iseqMetric = 1:2
            
            currSeqMetric = seqMetricMethods{iseqMetric};
            
            for iperiod = 1:3
                popBDscores.(periodNames{iperiod}).(currShuffleMethod).(currSeqMetric) = [popBDscores.(periodNames{iperiod}).(currShuffleMethod).(currSeqMetric); BDseqscore.(periodNames{iperiod}).data.(currShuffleMethod).(currSeqMetric).prctilescore(acceptedEvts.(periodNames{iperiod}))];  
            end
        end
    end

    
end


%% distributions


shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};
BDscoreMethods = {'weightedCorr'; 'replayScore'};
subfolder = '/home/kouroshmaboudi/Documents/HMM_project/GreatLakes_firstRun_Nov2020';

for ii = 1:numel(shuffleMethods)
    for jj = 1:numel(BDscoreMethods)

        shuffleMethod = shuffleMethods{ii};
        BDscoreMethod = BDscoreMethods{jj};
        
        
        figure;
        set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


        % % PRE
        subplot(3,3,[4 7])
        plotseqscoredist(popBDscores.PRE.(shuffleMethod).(BDscoreMethod), [], [], [], 0, sprintf('%s\n  (n=%d)', 'PRE', length(popBDscores.PRE.(shuffleMethod).(BDscoreMethod))))
        ylabel('ratio of PBEs', 'fontsize', 10)

        % % RUN
        subplot(3,3,[5 8])
        plotseqscoredist(popBDscores.RUN.(shuffleMethod).(BDscoreMethod), [], [], [], 0, sprintf('%s\n  (n=%d)', 'RUN', length(popBDscores.RUN.(shuffleMethod).(BDscoreMethod))))

        % % POST
        subplot(3,3,[6 9])
        plotseqscoredist(popBDscores.POST.(shuffleMethod).(BDscoreMethod), [], [], [], 0, sprintf('%s\n  (n=%d)', 'POST', length(popBDscores.POST.(shuffleMethod).(BDscoreMethod)))) 

        
        % legend
        legendSub = subplot(3,3,3);

        hold on

%         p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
        p_pts = plot(1, nan, 'color', 'r', 'linestyle', ':', 'linewidth', 2);
        p     = plot(1, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 3);

        h_legend = legend([p, p_pts],{'data', 'diagonal'}, 'Location', 'South');


        set(h_legend, 'fontsize', 8)
        legend boxoff 

        set(legendSub, 'Visible', 'off');
        
        
        if strcmp(BDscoreMethod, 'weightedCorr')
            namePart1 = 'weightedCorr';
        else
            namePart1 = 'radon transform';
        end
        
        
        filename = fullfile(subfolder, [BDscoreMethod '-' shuffleMethod]);
        
        print(gcf, filename, '-dpdf')
        


    end
end

close all

%%

subfolder = '/home/kouroshmaboudi/Documents/HMM_project/GreatLakes_firstRun_Nov2020';
plotScatterHistogram_RSvsWC_V2(popBDscores, subfolder)
plotScatterHistogram_TSvsPF_V2(popBDscores, subfolder)


%% a table of corr and p-values



sessionNames = fieldnames(popBDscores_bySession);
comparisons  = {'UI-RT vs UI-WC'; 'TS-RT vs TS-WC'; 'UI-RT vs TS-RT'; 'UI-WC vs TS-WC'; 'UI-RT vs TS-WC'; 'TS-RT vs UI-WC'};



dat = cell(numel(sessionNames), numel(comparisons));


for ii = 1:numel(sessionNames)
    
    
    % UI-RT vs UI-WC
    value1 = [popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.replayScore.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.replayScore.prctilescore];
    value2 = [popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.weightedCorr.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.weightedCorr.prctilescore];
    
    [corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
    if pval < 0.001
       dat{ii, 1} = sprintf('%.2f (p<0.001)', corrCoeff);
    else
       dat{ii, 1} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
    end
    
    
    % TS-RT vs TS-WC
    value1 = [popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.replayScore.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.replayScore.prctilescore];
    value2 = [popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.weightedCorr.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.weightedCorr.prctilescore];
    
    [corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
    if pval < 0.001
       dat{ii, 2} = sprintf('%.2f (p<0.001)', corrCoeff);
    else
       dat{ii, 2} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
    end
    
    
    % UI-RT vs TS-RT
    value1 = [popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.replayScore.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.replayScore.prctilescore];
    value2 = [popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.replayScore.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.replayScore.prctilescore];
    
    [corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
    if pval < 0.001
       dat{ii, 3} = sprintf('%.2f (p<0.001)', corrCoeff);
    else
       dat{ii, 3} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
    end
    
    
    % UI-WC vs TS-WC
    value1 = [popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.weightedCorr.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.weightedCorr.prctilescore];
    value2 = [popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.weightedCorr.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.weightedCorr.prctilescore];
    
    [corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
    if pval < 0.001
       dat{ii, 4} = sprintf('%.2f (p<0.001)', corrCoeff);
    else
       dat{ii, 4} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
    end
    
    
    
    % UI-RT vs TS-WC
    value1 = [popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.replayScore.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.replayScore.prctilescore];
    value2 = [popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.weightedCorr.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.weightedCorr.prctilescore];

    [corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
    if pval < 0.001
       dat{ii, 5} = sprintf('%.2f (p<0.001)', corrCoeff);
    else
       dat{ii, 5} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
    end
    
    
    % TS-RT vs UI-WC
    value1 = [popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.replayScore.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.wPBEtimeswap.replayScore.prctilescore];
    value2 = [popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.weightedCorr.prctilescore; popBDscores_bySession.(sessionNames{ii}).POST.data.unitIDshuffle.weightedCorr.prctilescore];
    
    [corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
    if pval < 0.001
       dat{ii, 6} = sprintf('%.2f (p<0.001)', corrCoeff);
    else
       dat{ii, 6} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
    end
    
end

% pooled result


% UI-RT vs UI-WC
value1 = [popBDscores.POST.unitIDshuffle.replayScore; popBDscores.POST.unitIDshuffle.replayScore];
value2 = [popBDscores.POST.unitIDshuffle.weightedCorr; popBDscores.POST.unitIDshuffle.weightedCorr];

[corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
if pval < 0.001
   dat{16, 1} = sprintf('%.2f (p<0.001)', corrCoeff);
else
   dat{16, 1} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
end

% TS-RT vs TS-WC
value1 = [popBDscores.POST.wPBEtimeswap.replayScore; popBDscores.POST.wPBEtimeswap.replayScore];
value2 = [popBDscores.POST.wPBEtimeswap.weightedCorr; popBDscores.POST.wPBEtimeswap.weightedCorr];

[corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
if pval < 0.001
   dat{16, 2} = sprintf('%.2f (p<0.001)', corrCoeff);
else
   dat{16, 2} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
end


% UI-RT vs TS-RT
value1 = [popBDscores.POST.unitIDshuffle.replayScore; popBDscores.POST.unitIDshuffle.replayScore];
value2 = [popBDscores.POST.wPBEtimeswap.replayScore; popBDscores.POST.wPBEtimeswap.replayScore];

[corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
if pval < 0.001
   dat{16, 3} = sprintf('%.2f (p<0.001)', corrCoeff);
else
   dat{16, 3} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
end


% UI-WC vs TS-WC
value1 = [popBDscores.POST.unitIDshuffle.weightedCorr; popBDscores.POST.unitIDshuffle.weightedCorr];
value2 = [popBDscores.POST.wPBEtimeswap.weightedCorr; popBDscores.POST.wPBEtimeswap.weightedCorr];

[corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
if pval < 0.001
   dat{16, 4} = sprintf('%.2f (p<0.001)', corrCoeff);
else
   dat{16, 4} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
end


% UI-RT vs TS-WC
value1 = [popBDscores.POST.unitIDshuffle.replayScore; popBDscores.POST.unitIDshuffle.replayScore];
value2 = [popBDscores.POST.wPBEtimeswap.weightedCorr; popBDscores.POST.wPBEtimeswap.weightedCorr];

[corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
if pval < 0.001
   dat{16, 5} = sprintf('%.2f (p<0.001)', corrCoeff);
else
   dat{16, 5} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
end


% TS-RT vs UI-WC
value1 = [popBDscores.POST.wPBEtimeswap.replayScore; popBDscores.POST.wPBEtimeswap.replayScore];
value2 = [popBDscores.POST.unitIDshuffle.weightedCorr; popBDscores.POST.unitIDshuffle.weightedCorr];

[corrCoeff, pval] = corr(value1, value2, 'type', 'spearman');
    
if pval < 0.001
   dat{16, 6} = sprintf('%.2f (p<0.001)', corrCoeff);
else
   dat{16, 6} = sprintf('%.2f (p=%.3f)', corrCoeff, pval);
end

sessionNames{16} = 'pooled';

figure;
uitable('columnname', comparisons, 'rowname', sessionNames, 'data', dat, 'columnwidth', {200}, 'position', [20 20 500 500]);




%% 

% plot all significance matrices corresponding to different jump distance
% criteria

shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};

sigMatrices = fieldnames(sigMatrix.RUN.data.(shuffleMethods{1})); % the different significance matrix 
jumpDistMatNames = {'max jump distance';'normalized max jump distance';'median jump distance'};

for ii = 1:numel(shuffleMethods)
    
    shuffleMethod = shuffleMethods{ii};
    
    for jj = 1: 3
        
        currSigMatrix = sigMatrices{jj};
        
       
        figure;
        set(gcf, 'Units', 'centimeters', 'position', [0 0 20 20])


        % PRE
        subplot(3,2,1);  plotsigMatrix(sigMatrix.PRE.data.(shuffleMethod).(currSigMatrix), {fileinfo.name;shuffleMethod;'';'Data'}, 0); xlabel(jumpDistMatNames{jj}); ylabel({'PRE', '','weighted correlation'});
        subplot(3,2,2);  plotsigMatrix(sigMatrix.PRE.p.(shuffleMethod).(currSigMatrix), 'pooled time swap', 0); xlabel(jumpDistMatNames{jj}); ylabel('weighted correlation'); 


        % RUN
        subplot(3,2,3);  plotsigMatrix(sigMatrix.RUN.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel({'RUN', '','weighted correlation'});
        subplot(3,2,4);  plotsigMatrix(sigMatrix.RUN.pts.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel('weighted correlation');


        % POST
        subplot(3,2,5);  plotsigMatrix(sigMatrix.POST.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel({'POST', '','weighted correlation'});
        subplot(3,2,6); plotsigMatrix(sigMatrix.POST.pts.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel('weighted correlation'); 


    
        filename = fullfile(directory, ['Bayesian__' currSigMatrix '__' shuffleMethod]);
    
        savepdf(gcf, filename, '-dpng')
%         saveas(gcf, filename,  'epsc')
    end
    
    jj = 4;
    
    currSigMatrix = sigMatrices{jj};
    figure;
    set(gcf, 'Units', 'centimeters', 'position', [0 0 20 20])


    % PRE
    subplot(3,2,1);  plotsigMatrix_rs(sigMatrix.PRE.data.(shuffleMethod).(currSigMatrix), {fileinfo.name; shuffleMethod;'';'Data'}, 0); xlabel('covered track length'); ylabel({'PRE', '','replay score'});
    subplot(3,2,2);  plotsigMatrix_rs(sigMatrix.PRE.pts.(shuffleMethod).(currSigMatrix), 'poisson', 0); xlabel('covered track length'); ylabel('replay score'); 


    % RUN
    subplot(3,2,3);  plotsigMatrix_rs(sigMatrix.RUN.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel({'RUN', '','replay score'});
    subplot(3,2,4);  plotsigMatrix_rs(sigMatrix.RUN.pts.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel('replay score');


    % POST
    subplot(3,2,5);  plotsigMatrix_rs(sigMatrix.POST.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel({'POST', '','replay score'});
    subplot(3,2,6); plotsigMatrix_rs(sigMatrix.POST.pts.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel('replay score'); 



    filename = fullfile(directory, ['Bayesian__' currSigMatrix '__' shuffleMethod]);

    savepdf(gcf, filename, '-dpng')
%     saveas(gcf, filename,  'epsc')


end







