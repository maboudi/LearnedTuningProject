




nUnits = size(spatialTunings_merge, 1);

shifts = -floor(nPosBins/2):1:floor(nPosBins/2);
numShifts = numel(shifts);


corrCoeff_zeroShift = zeros(nUnits, 1);
KLdiv_zeroShift     = zeros(nUnits, 1);

bb = spatialTunings_merge(PFunit, :);
bb = bb/sum(bb);

for iUnit = 1:nUnits
    
    LTunit = iUnit;
    
    aa = currLT(LTunit, :);
    aa = aa/sum(aa);

    
    
    corrCoeff = zeros(numShifts, 1);
    KLdiv     = zeros(numShifts, 1);
    
    for ii = 1:numel(shifts)

        currShift = shifts(ii);
        aa2 = circshift(aa, currShift);

        corrCoeff(ii) = corr(aa2', bb', 'type', 'pearson');
        KLdiv(ii)     = nansum(bb .* (log2(bb) - log2(aa2)));

    end
    corrCoeff_norm = (corrCoeff - mean(corrCoeff))/std(corrCoeff);
    KLdiv_norm     = (KLdiv - mean(KLdiv))/std(KLdiv);

    corrCoeff_zeroShift(iUnit) = mean(corrCoeff_norm(shifts > -5 & shifts < 5));
    KLdiv_zeroShift(iUnit)     = mean(KLdiv_norm(shifts > -5 & shifts < 5));

    
    
end

%%
close all

PFunit = 62;
LTunit = 62;

figure;

customPlot(currLT(LTunit,:), spatialTunings_merge(PFunit, :), PFunit, LTunit, iEpoch)





function customPlot(learnedTuning, spatialTuning, PFunit, LTunit, iEpoch)



epochNames = {'PRE'; 'RUN'; 'POST'};

subplot(2,1,1)



aa = learnedTuning;
bb = spatialTuning;

aa = aa/sum(aa);
bb = bb/sum(bb);

KLbeforeSum = bb .* (log2(bb) - log2(aa));
KLdiv = nansum(KLbeforeSum); 


nPosBins = numel(aa);

hold on
yyaxis right
plot((1:nPosBins)*2, aa, 'DisplayName', 'learned tuning');
ylabel('learend tuning', 'fontsize', 8)

yyaxis left
plot((1:nPosBins)*2, bb, 'DisplayName', 'place field')
ylabel('place field', 'fontsize', 8)
xlabel('position (cm)', 'fontsize', 8)


title(sprintf('place field of unit %d and %s Learned tuning of unit %d', PFunit, epochNames{iEpoch}, LTunit), 'fontsize', 8)

grid on

xl = xlim;
yl = ylim;

corrCoeff = corr(aa', bb', 'type', 'pearson');

text(median(xl), yl(2)-0.05*range(yl), sprintf('Pearson corr. = %.2f', corrCoeff), 'fontsize', 8)
text(median(xl), yl(2)-0.15*range(yl), sprintf('KLdiv. = %.2f', KLdiv), 'fontsize', 8)


legend('location', 'southeast', 'box', 'off')



%%%

subplot(2,1,2)

shifts = -floor(nPosBins/2):1:floor(nPosBins/2);

corrCoeff = zeros(numel(shifts), 1);
KLdiv     = zeros(numel(shifts), 1);

for ii = 1:numel(shifts)
    
    currShift = shifts(ii);
    aa2 = circshift(aa, currShift);
   
    corrCoeff(ii) = corr(aa2', bb', 'type', 'pearson');
    KLdiv(ii)     = nansum(bb .* (log2(bb) - log2(aa2)));
    
end

hold on

yyaxis left
plot(2*shifts, corrCoeff)
ylabel('Pearson correlation', 'fontsize', 8)

yyaxis right
plot(2*shifts, KLdiv)
ylabel('KL divergence', 'fontsize', 8)

xlabel('position shift', 'fontsize', 8)

grid on



corrCoeff_norm = (corrCoeff - mean(corrCoeff))/std(corrCoeff);
KLdiv_norm     = (KLdiv - mean(KLdiv))/std(KLdiv);


corrCoeff_zeroShift = mean(corrCoeff_norm(shifts > -5 & shifts < 5));
KLdiv_zeroShift     = mean(KLdiv_norm(shifts > -5 & shifts < 5));

yl = ylim;

text(0, yl(2)-0.05*range(yl), sprintf('corr(norm) = %.2f', corrCoeff_zeroShift), 'fontsize', 8)
text(0, yl(2)-0.15*range(yl), sprintf('KLdiv(norm) = %.2f', KLdiv_zeroShift), 'fontsize', 8)



[maxCorrCoeff, maxShiftIdx_c] = max(corrCoeff_norm);
[maxKLdiv, maxShiftIdx_k] = min(KLdiv_norm);


if (shifts(maxShiftIdx_c) < -5 || shifts(maxShiftIdx_c) > 5) && maxCorrCoeff > 2
    yl = ylim;
    text(2*shifts(maxShiftIdx_c), median(yl), sprintf('corr(norm) = %.2f, shift = %d', maxCorrCoeff, 2*shifts(maxShiftIdx_c)), 'fontsize', 8)
end

if (shifts(maxShiftIdx_k) < -5 || shifts(maxShiftIdx_k) > 5) && maxKLdiv < -2
    yl = ylim;
    text(2*shifts(maxShiftIdx_k), median(yl)-0.1*range(yl), sprintf('KLdiv(norm) = %.2f, shift = %d', maxKLdiv, 2*shifts(maxShiftIdx_k)), 'fontsize', 8)

end

end
