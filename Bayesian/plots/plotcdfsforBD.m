clc;
close all


subfolder = fullfile(mainDir, 'Bayesian'); 



%% PRE

% data

load(fullfile(subfolder, 'PRE', 'data', 'BDresults_wc_jd.mat'))

actualPRE = BDseqscore.zscore;

% actualPRE = BDprctile(:, 1, :); % just the time swap
% actualPRE = squeeze(actualPRE);
actualPRE = max(actualPRE, [], 2);

% shuffle

load(fullfile(subfolder, 'PRE', 'withinPBETimeSwap', 'BDresults_wc_jd.mat'))

shufflePRE = BDseqscore.zscore;

% shufflePRE = BDprctile(:, 1, :);
% shufflePRE = squeeze(shufflePRE);
shufflePRE = max(shufflePRE, [], 2);



%% RUN


% data

load(fullfile(subfolder, 'RUN', 'data', 'BDresults_wc_jd.mat'))

actualRUN = BDseqscore.zscore;

% actualPRE = BDprctile(:, 1, :); % just the time swap
% actualPRE = squeeze(actualPRE);
actualRUN = max(actualRUN, [], 2);

% shuffle

load(fullfile(subfolder, 'RUN', 'withinPBETimeSwap', 'BDresults_wc_jd.mat'))

shuffleRUN = BDseqscore.zscore;

% shuffleRUN = BDprctile(:, 1, :);
% shuffleRUN = squeeze(shuffleRUN);
shuffleRUN = max(shuffleRUN, [], 2);



%% POST


% data

load(fullfile(subfolder, 'POST', 'data', 'BDresults_wc_jd.mat'))

actualPOST = BDseqscore.zscore;

% actualPOST = BDprctile(:, 1, :); % just the time swap
% actualPOST = squeeze(actualPOST);
actualPOST = max(actualPOST, [], 2);

% shuffle

load(fullfile(subfolder, 'POST', 'withinPBETimeSwap', 'BDresults_wc_jd.mat'))

shufflePOST = BDseqscore.zscore;

% shufflePOST = BDprctile(:, 1, :);
% shufflePOST = squeeze(shufflePOST);
shufflePOST = max(shufflePOST, [], 2);



%% plot


figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


% Cross-validations 

subplot(1,3,1)
congruencePlot2(actualPRE/100, shufflePRE/100, 1, sprintf('%s\n  (n=%d)', 'PRE', length(actualPRE)))

subplot(1,3,2)
congruencePlot2(actualRUN/100, shuffleRUN/100, 0, sprintf('%s\n  (n=%d)', 'RUN', length(actualRUN)))

subplot(1,3,3)
congruencePlot2(actualPOST/100, shufflePOST/100, 0, sprintf('%s\n  (n=%d)', 'POST', length(actualPOST)))


filename = fullfile(subfolder, 'BDcdfs');

savepdf(gcf, filename, '-dsvg')
savepdf(gcf, filename, '-dpng')



