
% 
% tuningCurves = spatialTunings_LR2;
% 
% [peaks, peakPositions] = max(tuningCurves, [], 2);
% [~, sortIdx] = sort(peakPositions, 'ascend');
% 
% 
% tuningCurves = tuningCurves(sortIdx, :);
% 
% tuningCurves2 = tuningCurves ./ repmat(max(tuningCurves, [], 2), [1 size(tuningCurves, 2)]);
% 
% figure; imagesc(tuningCurves2); colormap('jet'); 


posteriorsPRE  = posteriorProbMatrix.PRE.data(:, 1); conctPosteriorsPRE = cell2mat(posteriorsPRE'); sharpnessPRE = max(conctPosteriorsPRE,[], 1);
posteriorsRUN  = posteriorProbMatrix.RUN.data(:, 1); conctPosteriorsRUN = cell2mat(posteriorsRUN'); sharpnessRUN = max(conctPosteriorsRUN,[], 1);
posteriorsPOST = posteriorProbMatrix.POST.data(:, 1); conctPosteriorsPOST = cell2mat(posteriorsPOST'); sharpnessPOST = max(conctPosteriorsPOST,[], 1);


concatSharpness = [sharpnessPRE sharpnessRUN sharpnessPOST]';


conctFiringsPRE  = cell2mat(PREbinnedPBEs.data(:, 2)'); conctFiringsPRE(find(conctFiringsPRE)) = 1; nUnitsFiringPRE = sum(conctFiringsPRE, 1);
conctFiringsRUN  = cell2mat(RUNbinnedPBEs.data(:, 2)'); conctFiringsRUN(find(conctFiringsRUN)) = 1; nUnitsFiringRUN = sum(conctFiringsRUN, 1);
conctFiringsPOST = cell2mat(POSTbinnedPBEs.data(:, 2)'); conctFiringsPOST(find(conctFiringsPOST)) = 1; nUnitsFiringPOST = sum(conctFiringsPOST, 1);

concatnUnitsFiring = [nUnitsFiringPRE nUnitsFiringRUN nUnitsFiringPOST]';

idx = find(concatnUnitsFiring);

concatSharpness = concatSharpness(idx);
concatnUnitsFiring = concatnUnitsFiring(idx);

figure; 

subplot(3, 2, [1 3])

randSelect =randperm(length(concatSharpness), 5000);
plot(concatnUnitsFiring(randSelect), concatSharpness(randSelect), '.', 'markersize', 5)

rho = corr(concatnUnitsFiring', concatSharpness');


c = polyfit(concatnUnitsFiring, concatSharpness,1);
sharpness_est = polyval(c, concatnUnitsFiring);


currPos = get(gca, 'position');

set(gca, 'position', [currPos(1) currPos(2)+0.1 currPos(3) currPos(4)-0.2])

hold on


plot(concatnUnitsFiring(randSelect)', sharpness_est(randSelect)', 'linewidth', 2, 'color', 'r')

xlabel('number of firing units', 'fontsize', 10)
ylabel({'peak posterior pr'; '(sharpness)'}, 'fontsize', 10)

ylim([0 0.7])
xlim([0 13])

text(2, 0.6, sprintf('rho=%.2f', rho), 'fontsize', 10, 'color', 'r')


subplot(3, 2, 5)

bar([median(sharpnessPRE) median(sharpnessRUN) median(sharpnessPOST)])
set(gca, 'xticklabel', {'PRE', 'RUN', 'POST'}, 'box', 'off')
ylabel('median sharpness', 'fontsize', 10)

hold on
errorbar([1 2 3],[median(sharpnessPRE) median(sharpnessRUN) median(sharpnessPOST)],[prctile(sharpnessPRE, 25) prctile(sharpnessRUN, 25) prctile(sharpnessPOST, 25)], [prctile(sharpnessPRE, 75) prctile(sharpnessRUN, 75) prctile(sharpnessPOST, 75)])


ax1 = subplot(3,2, [2 4]);

bins = 0:0.1:max([sharpnessPRE sharpnessPOST]);

countPRE = hist(sharpnessPRE, bins); countPRE = countPRE/sum(countPRE);
countRUN = hist(sharpnessRUN, bins); countRUN = countRUN/sum(countRUN);
countPOST = hist(sharpnessPOST, bins); countPOST = countPOST/sum(countPOST);



bar(bins, countPRE, 'faceColor', 'b', 'edgeColor', 'none', 'facealpha', 0.5)
hold on
bar(bins, countPOST, 'faceColor', 'r', 'edgeColor', 'none', 'facealpha', 0.5)
hold on
bar(bins, countRUN, 'faceColor', 'none', 'edgeColor', 'k')


legend('PRE', 'POST', 'RUN')


currPos = get(gca, 'position');
set(gca, 'position', [currPos(1) currPos(2)+0.1 currPos(3) currPos(4)-0.2])

xlabel('sharpness', 'fontsize', 10)
ylabel('probability')

% ylim([0.1 0.13])

subplot(3,2,6)

bar([median(nUnitsFiringPRE) median(nUnitsFiringRUN) median(nUnitsFiringPOST)])
set(gca, 'xticklabel', {'PRE', 'RUN', 'POST'}, 'box', 'off')
ylabel('median number of firing units', 'fontsize', 10)


hold on
errorbar([1 2 3],[median(nUnitsFiringPRE) median(nUnitsFiringRUN) median(nUnitsFiringPOST)],[prctile(nUnitsFiringPRE, 25) prctile(nUnitsFiringRUN, 25) prctile(nUnitsFiringPOST, 25)], [prctile(nUnitsFiringPRE, 75) prctile(nUnitsFiringRUN, 75) prctile(nUnitsFiringPOST, 75)])



dim = [0 0.9 0.2 .1];
textOnTop = {'Sharpness of the posterior probability versus number of units firing within individual time bins'; ...
    'Session: Achilles_10252013'};
    
annotation('textbox',dim,'String',textOnTop, 'FitBoxToText','on', 'Interpreter', 'none');


print(gcf, 'sharpnessVSnumFiringUnits', '-dpdf', '-r0')

