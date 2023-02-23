clc;
close all

directory = fullfile(mainDir, 'Bayesian'); 


%% PRE


% data
actualPRE  = loadfile('PRE', 'data_twoShuffleMethods', directory); 

% poisson
poissonPRE = loadfile('PRE', 'poisson', directory); 

% pooled time swap
ptsPRE     = loadfile('PRE', 'pooledTimeSwap', directory);

% within PBE TimeSwap
wtsPRE     = loadfile('PRE', 'withinPBETimeSwap', directory);

% 
% 
% 
% %% RUN
% 
% 
% % data
% actualRUN  = loadfile('RUN', 'data_twoShuffleMethods', directory); 
% 
% % poisson
% poissonRUN = loadfile('RUN', 'poisson', directory); 
% 
% % pooled time swap
% ptsRUN     = loadfile('RUN', 'pooledTimeSwap', directory);
% 
% % within PBE TimeSwap
% wtsRUN     = loadfile('RUN', 'withinPBETimeSwap', directory);
% 
% 
% 
% %% POST
% 
% 
% % data
% actualPOST  = loadfile('RUN', 'data_twoShuffleMethods', directory); 
% 
% % poisson
% poissonPOST = loadfile('RUN', 'poisson', directory); 
% 
% % pooled time swap
% ptsPOST     = loadfile('RUN', 'pooledTimeSwap', directory);
% 
% % within PBE TimeSwap
% wtsPOST     = loadfile('RUN', 'withinPBETimeSwap', directory);



%% plot


figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


% Cross-validations 

subplot(1,3,1)
okLegend = 1;
congruencePlot(actualPRE, poissonPRE, wtsPRE, ptsPRE, okLegend, sprintf('%s\n  (n=%d)', 'PRE', length(actualPRE)))


% subplot(1,3,2)
% congruencePlot2(actualRUN, shuffleRUN, 0, sprintf('%s\n  (n=%d)', 'RUN', length(actualRUN)))
% 
% subplot(1,3,3)
% congruencePlot2(actualPOST, shufflePOST, 0, sprintf('%s\n  (n=%d)', 'POST', length(actualPOST)))


filename = fullfile(directory, 'BDcdfs_2');

savepdf(gcf, filename, '-dsvg')
savepdf(gcf, filename, '-dpng')





%% functions

function data = loadfile(period, dataType, directory)

load(fullfile(directory, period, dataType, 'BDresults.mat'))

data = BDprctile(:, 1, :); % just the time swap
data = squeeze(data);
data = max(data, [], 2);

end


