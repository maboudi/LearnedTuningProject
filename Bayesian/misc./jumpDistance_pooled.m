clear
clc

currDir = '~/Documents/HMM_project/GreatLakes_March2021_addingMoreShuffles_CircularPF_CircularPosterior';
cd(currDir)
rr = dir;

% periodNames = {'PRE'; 'RUN'; 'POST'};
periodNames = {'RUN'};

shuffleMethodNames = {'wPBEtimeswap';'unitIDshuffle'; 'circularPFshuffle'; 'circularPosteriorShuffle'};
seqMetricNames = {'weightedCorr', 'RadonTF'};


jumpDist_pooled      = struct('PRE', [], 'RUN', [], 'POST', []);
jumpDist_null_pooled = struct('PRE', [], 'RUN', [], 'POST', []);

weightedCorr_pooled = struct('PRE', [], 'RUN', [], 'POST', []);
weightedCorr_null_pooled = struct('PRE', [], 'RUN', [], 'POST', []);

RadonTF_pooled = struct('PRE', [], 'RUN', [], 'POST', []);
RadonTF_null_pooled = struct('PRE', [], 'RUN', [], 'POST', []);

coveredLen_pooled = struct('PRE', [], 'RUN', [], 'POST', []);
coveredLen_null_pooled = struct('PRE', [], 'RUN', [], 'POST', []);



for iperiod = 1:numel(periodNames)
    for ishuffle = 1:4
        for iseqMetric = 1:2

            BDseqscore_pooled.(periodNames{iperiod}).(shuffleMethodNames{ishuffle}).(seqMetricNames{iseqMetric}) = [];

        end
    end
end



for sessionNumber = 3

    
    
    
    sessionName = rr(sessionNumber+2).name

    currSessionFolder = fullfile(currDir, sessionName);

    cd(currSessionFolder)

    cd('poisson')
    
    load('Bayesian.mat', 'BDseqscore', 'jumpDist', 'jumpDist_null', 'weightedCorr', 'weightedCorr_null', 'RadonTF', 'RadonTF_null', 'coveredLen', 'coveredLen_null', ...
        'secondaryPBEs', 'behavior', 'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs')
    
    
    datasets.PRE = PREbinnedPBEs;
    datasets.RUN = RUNbinnedPBEs;
    datasets.POST = POSTbinnedPBEs;


    PBELen = struct('PRE', [], 'RUN', [], 'POST', []);

    for pp = 1:numel(periodNames)

        currentPBEs = datasets.(periodNames{pp});
        nPBEs = size(currentPBEs, 1);

        PBELen.(periodNames{pp}) = zeros(nPBEs, 1);
        nFiringUnits.(periodNames{pp}) = zeros(nPBEs, 1);

        for pbe = 1:nPBEs

            PBELen.(periodNames{pp})(pbe) = size(currentPBEs{pbe, 2}, 2);

            nFiringUnits.(periodNames{pp})(pbe) = numel(find(sum(currentPBEs{pbe, 2}, 2)));

        end

    end
    
    
    try
        secondaryPBEs_1.PRE = secondaryPBEs(secondaryPBEs(:, 1) > behavior.PREEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.PREEpoch(1, 2), :);
        secondaryPBEs_1.RUN = secondaryPBEs(secondaryPBEs(:, 1) > behavior.MazeEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.MazeEpoch(1, 2) & secondaryPBEs(:, 5) == 1, :);
        secondaryPBEs_1.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.POSTEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.POSTEpoch(1, 2), :);
    catch
        secondaryPBEs_1.PRE = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(1, 1) & secondaryPBEs(:, 2) < behavior.time(1, 2), :);
        secondaryPBEs_1.RUN = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(2, 1) & secondaryPBEs(:, 2) < behavior.time(2, 2) & secondaryPBEs(:, 5) == 1, :);
        secondaryPBEs_1.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.time(3, 1) & secondaryPBEs(:, 2) < behavior.time(3, 2), :);
    end

%     acceptedEvts.PRE = find(sum(secondaryPBEs_1.PRE(:, [6]), 2) == 1 & PBELen.PRE >= 4 & nFiringUnits.PRE >= 5);
    acceptedEvts.RUN = find(secondaryPBEs_1.RUN(:, 5) == 1 & PBELen.RUN >= 4 & nFiringUnits.RUN >= 5);
%     acceptedEvts.POST = find(sum(secondaryPBEs_1.POST(:, [6]), 2) == 1 & PBELen.POST >= 4 & nFiringUnits.POST >= 5);

    
    for ii = 1:numel(periodNames)
        
        currPeriod = periodNames{ii};
        


        currJumpDistances = [jumpDist.(currPeriod).p(acceptedEvts.(currPeriod))]';
        jumpDist_pooled.(currPeriod) = [jumpDist_pooled.(currPeriod); currJumpDistances];

        currJumpDistances_null = [jumpDist_null.(currPeriod).p(acceptedEvts.(currPeriod))]';
        jumpDist_null_pooled.(currPeriod) = [jumpDist_null_pooled.(currPeriod); currJumpDistances_null];


        weightedCorr_pooled.(currPeriod) = [weightedCorr_pooled.(currPeriod); [weightedCorr.(currPeriod).p(acceptedEvts.(currPeriod))]];
        weightedCorr_null_pooled.(currPeriod) = [weightedCorr_null_pooled.(currPeriod); [weightedCorr_null.(currPeriod).p(acceptedEvts.(currPeriod), :, :)]];


        RadonTF_pooled.(currPeriod) = [RadonTF_pooled.(currPeriod); [RadonTF.(currPeriod).p(acceptedEvts.(currPeriod))]];
        RadonTF_null_pooled.(currPeriod) = [RadonTF_null_pooled.(currPeriod); [RadonTF_null.(currPeriod).p(acceptedEvts.(currPeriod), :, :)]];


        coveredLen_pooled.(currPeriod) = [coveredLen_pooled.(currPeriod); [coveredLen.(currPeriod).p(acceptedEvts.(currPeriod))]];
        coveredLen_null_pooled.(currPeriod) = [coveredLen_null_pooled.(currPeriod); [coveredLen_null.(currPeriod).p(acceptedEvts.(currPeriod), :, :)]];
        
        
        for ishuffle = 1:4
            for iseqMetric = 1:2
                
                BDseqscore_pooled.(currPeriod).(shuffleMethodNames{ishuffle}).(seqMetricNames{iseqMetric}) = [BDseqscore_pooled.(currPeriod).(shuffleMethodNames{ishuffle}).(seqMetricNames{iseqMetric}), ...
                                                                                                                [BDseqscore.(currPeriod).p.(shuffleMethodNames{ishuffle}).(seqMetricNames{iseqMetric}).prctilescore(acceptedEvts.(currPeriod))]'];
                
            end
        end
        
        
        
        
    end


end




%%


nShuffles = 500;



for ip = 1:numel(periodNames)
    
    currPeriod = periodNames{ip};

    nEvents = length(jumpDist_pooled.(currPeriod));


    for shuffleMethod = 1:4


        %%% absolute weighted correlation and jump distance


        maxJump     = zeros(nEvents, 1);
        normMaxJump = zeros(nEvents, 1);
        medianJump  = zeros(nEvents, 1);

        for pbe = 1:nEvents    
            maxJump(pbe)     = jumpDist_pooled.(currPeriod)(pbe).maxJump;
            normMaxJump(pbe) = jumpDist_pooled.(currPeriod)(pbe).normMaxJump;
            medianJump(pbe)  = jumpDist_pooled.(currPeriod)(pbe).medianJump;
        end

        maxJump_surr     = zeros(nEvents, nShuffles);
        normMaxJump_surr = zeros(nEvents, nShuffles);
        medianJump_Surr  = zeros(nEvents, nShuffles);

        for pbe = 1: nEvents 

          maxJump_surr(pbe, :)     = jumpDist_null_pooled.(currPeriod)(pbe).maxJump(:, shuffleMethod);
          normMaxJump_surr(pbe, :) = jumpDist_null_pooled.(currPeriod)(pbe).normMaxJump(:, shuffleMethod);
          medianJump_Surr(pbe, :)  = jumpDist_null_pooled.(currPeriod)(pbe).medianJump(:, shuffleMethod);

        end

        % calculating the matrices of PBE counts passing pairs of threshold on the two
        % metrics

        wcthresholds = 0:0.1:0.9;
        jdthresholds = 0.1:0.1:1;

        nwcthreshs   = length(wcthresholds);
        njdthreshs   = length(jdthresholds);



        % DATA

        % using max jump distance
        twoDMatrix1 = calRScoreJdistMatrix(weightedCorr_pooled.(currPeriod), maxJump, wcthresholds, jdthresholds);

        % using normalized max jump distance
        twoDMatrix2 = calRScoreJdistMatrix(weightedCorr_pooled.(currPeriod), normMaxJump, wcthresholds, jdthresholds);

        % using median jump distance
        twoDMatrix3 = calRScoreJdistMatrix(weightedCorr_pooled.(currPeriod), medianJump, wcthresholds, jdthresholds);



        % SURROGATES

        twoDMatrix_surr1 = zeros(nwcthreshs, njdthreshs, nShuffles); % max jump distance
        twoDMatrix_surr2 = zeros(nwcthreshs, njdthreshs, nShuffles); % normalized max jump distance
        twoDMatrix_surr3 = zeros(nwcthreshs, njdthreshs, nShuffles); % median Jump distance

        for sn = 1: nShuffles

            % max jump distance
            twoDMatrix_surr1(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null_pooled.(currPeriod)(:, sn, shuffleMethod), maxJump_surr(:, sn), wcthresholds, jdthresholds);

            % normalized max jump distance
            twoDMatrix_surr2(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null_pooled.(currPeriod)(:, sn, shuffleMethod), normMaxJump_surr(:, sn), wcthresholds, jdthresholds);

            % median Jump distance
            twoDMatrix_surr3(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null_pooled.(currPeriod)(:, sn, shuffleMethod), medianJump_Surr(:, sn), wcthresholds, jdthresholds);

        end


        % SIGNIFICANCE MATRIX

        sigMatrix.(currPeriod).(shuffleMethodNames{shuffleMethod}) = struct('wc_maxJump', zeros(nwcthreshs, njdthreshs), 'wc_normMaxJump', zeros(nwcthreshs, njdthreshs), 'wc_medianJump', zeros(nwcthreshs, njdthreshs), 'rt_slp', []);


        for ii = 1: nwcthreshs
           for jj = 1: njdthreshs

               nSurrGreater = length(find(twoDMatrix_surr1(ii, jj, :) >= twoDMatrix1(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(currPeriod).(shuffleMethodNames{shuffleMethod}).wc_maxJump(ii, jj) = nSurrGreater/nShuffles; % Monte Carlo P-value

               nSurrGreater = length(find(twoDMatrix_surr2(ii, jj, :) >= twoDMatrix2(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(currPeriod).(shuffleMethodNames{shuffleMethod}).wc_normMaxJump(ii, jj) = nSurrGreater/nShuffles;

               nSurrGreater = length(find(twoDMatrix_surr3(ii, jj, :) >= twoDMatrix3(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(currPeriod).(shuffleMethodNames{shuffleMethod}).wc_medianJump(ii, jj) = nSurrGreater/nShuffles;


           end
        end


        %% significance matrix2: replay score and coveredLen


        %%% significance matrix

        rsthresholds = 0:0.1:0.9;
        coveredLenthresholds = 0:0.1:0.9; 

        nrsthreshs         = length(rsthresholds);
        ncoveredLenthreshs = length(coveredLenthresholds);


        twoDMatrix = calRScoreJdistMatrix2(RadonTF_pooled.(currPeriod), coveredLen_pooled.(currPeriod), rsthresholds, coveredLenthresholds);


        % SURROGATES

        twoDMatrix_surr = zeros(nrsthreshs, ncoveredLenthreshs, nShuffles); % max jump distance

        for sn = 1: nShuffles   
            twoDMatrix_surr(:, :, sn) = calRScoreJdistMatrix2(RadonTF_null_pooled.(currPeriod)(:, sn, shuffleMethod), coveredLen_null_pooled.(currPeriod)(:, sn, shuffleMethod), rsthresholds, coveredLenthresholds);
        end


        sigMatrix.(currPeriod).(shuffleMethodNames{shuffleMethod}).rt_slp = zeros(nrsthreshs, ncoveredLenthreshs);
        for ii = 1: nrsthreshs
           for jj = 1: ncoveredLenthreshs

               nSurrGreater     = length(find(twoDMatrix_surr(ii, jj, :) >= twoDMatrix(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(currPeriod).(shuffleMethodNames{shuffleMethod}).rt_slp(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

           end
        end



    end

end



%% plot the 2-D significance matrices

figure; 

set(gcf, 'units', 'centimeter', 'position', [13 10 22 21])

sigMatrixNames = {'wc_maxJump';'wc_normMaxJump'; 'wc_medianJump';'rt_slp'};

imat = 3;
currMat = sigMatrixNames{imat};
superttl = {'weighted corr - maximum jump distance 2-factor pval matrix'; ...
    'weighted corr - normalized maximum jump distance 2-factor pval matrix'; ...
    'weighted corr - median jump distance 2-factor pval matrix'; ...
    'radon transfrom - track covered length 2-factor pval matrix'};

shuffleMethodsFullNames = {...
    'w/ PBE time swap', ...
    'unit id shuffle', ...
    'place field c. shuffle', ...
    'posterior dist. c. shuffle'};

for ism = 1:numel(shuffleMethodNames)
    
    currShuffleMethod = shuffleMethodNames{ism};
    
    for iperiod = 1:numel(periodNames)
        
        currPeriod = periodNames{iperiod};
        
        rr = subplot(4,3,(ism-1)*3+iperiod);
        
        
        if ism == 1 && iperiod == 3
            OKcolorbar = 1;
        else
            OKcolorbar = 0;
        end
        
        
        if imat <=3
            plotsigMatrix(sigMatrix.(currPeriod).(currShuffleMethod).(currMat), '', OKcolorbar)
        else
            plotsigMatrix_rs(sigMatrix.(currPeriod).(currShuffleMethod).(currMat), '', OKcolorbar)
        end
        
        if ism == 1
            title(currPeriod, 'fontsize', 14, 'fontweight', 'normal')
        end
        
        
        if ism == numel(shuffleMethodNames) && iperiod == 1
            
            switch imat
                case 1
                    xlabel('max jump distance')
                    ylabel('abs w. correlation')
                case 2
                    xlabel('norm. max jump distance')
                    ylabel('abs w. correlation')
                case 3
                    xlabel('median jump distance')
                    ylabel('abs w. correlation')
                case 4
                    xlabel('covered length')
                    ylabel('radon integral')
            end
                    
        end
        
        if iperiod == 1
            oldYLabel = rr.YLabel.String;            
            ylabel({shuffleMethodsFullNames{ism}; ''; oldYLabel})
        end
            
    end
end
suptitle(superttl{imat})



%%

seqScoreMethod = 'weightedCorr';
        
        
figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])

% % PRE

% subplot(3,3,[4 7])
% plotseqscoredist(...
%     BDseqscore_pooled.PRE.wPBEtimeswap.(seqScoreMethod), ... 
%     BDseqscore_pooled.PRE.unitIDshuffle.(seqScoreMethod), ...
%     BDseqscore_pooled.PRE.circularPFshuffle.(seqScoreMethod), ...
%     BDseqscore_pooled.PRE.circularPosteriorShuffle.(seqScoreMethod), ...
%     0, 'PRE')
% 
% 
% ylabel('ratio of PBEs', 'fontsize', 10)

% % RUN
subplot(3,3,[5 8])
plotseqscoredist(...
    BDseqscore_pooled.RUN.wPBEtimeswap.(seqScoreMethod), ... 
    BDseqscore_pooled.RUN.unitIDshuffle.(seqScoreMethod), ...
    BDseqscore_pooled.RUN.circularPFshuffle.(seqScoreMethod), ...
    BDseqscore_pooled.RUN.circularPosteriorShuffle.(seqScoreMethod), ...
    0, 'RUN')
% % % POST
% subplot(3,3,[6 9])
% plotseqscoredist(...
%     BDseqscore_pooled.POST.wPBEtimeswap.(seqScoreMethod), ... 
%     BDseqscore_pooled.POST.unitIDshuffle.(seqScoreMethod), ...
%     BDseqscore_pooled.POST.circularPFshuffle.(seqScoreMethod), ...
%     BDseqscore_pooled.POST.circularPosteriorShuffle.(seqScoreMethod), ...
%     0, 'POST')

% legend
legendSub = subplot(3,3,3);

hold on

p_ts        = plot(1, nan, 'color', [0.4940    0.1840    0.5560 0.5], 'linestyle', '-', 'linewidth', 2);
p_ui        = plot(1, nan, 'color', [0.9290    0.6940    0.1250 0.5], 'linestyle', '-', 'linewidth', 2);
p_pf        = plot(1, nan, 'color', [0.8500    0.3250    0.0980 0.5], 'linestyle', '-', 'linewidth', 2);
p_posterior = plot(1, nan, 'color', [0         0.4470    0.7410 0.5], 'linestyle', '-', 'linewidth', 2);


h_legend = legend([p_ts, p_ui, p_pf, p_posterior],{'time swap', 'unit id shuffle', 'place field c. shuffle', 'posterior c. shuffle'}, 'Location', 'South');


set(h_legend, 'fontsize', 8)
legend boxoff 

set(legendSub, 'Visible', 'off');


%%


function count = calRScoreJdistMatrix(weightedCorr, jumpDistance, wcthresholds, jdthresholds)

nwcthreshs = length(wcthresholds);
njdthreshs = length(jdthresholds);

count = zeros(nwcthreshs, njdthreshs);

for ii = 1:nwcthreshs
    
    wcthresh = wcthresholds(ii);
    
    for jj = 1: njdthreshs
        
        jdthresh = jdthresholds(jj);
        
        count(ii, jj) = length(find(abs(weightedCorr) >= wcthresh & jumpDistance <= jdthresh)); % number of PBEs meeting passing the two thresholds
        
    end
end

end


function count = calRScoreJdistMatrix2(RadonTF, coveredLen, rsthresholds, coveredLenthresholds)

nrsthreshs = length(rsthresholds);
ncoveredLenthreshs = length(coveredLenthresholds);

count = zeros(nrsthreshs, ncoveredLenthreshs);

for ii = 1:nrsthreshs
    
    rsthresh = rsthresholds(ii);
    
    for jj = 1: ncoveredLenthreshs
        
        coveredLenthresh = coveredLenthresholds(jj);
        
        count(ii, jj) = length(find(RadonTF >= rsthresh & coveredLen >= coveredLenthresh)); % number of PBEs meeting passing the two thresholds
        
    end
end


end



