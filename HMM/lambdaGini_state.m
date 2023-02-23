function [realGiniCoeff, unitCycleGiniCoeff,timeSwapGiniCoeff, poissGiniCoeff, identityGiniCoeff] = lambdaGini_state(lambda, lambda_c, lambda_ts, lambda_p, lambda_i, ii)

real = lambda(:,:,ii);
unitCycle = lambda_c(:,:,ii);
timeSwap = lambda_ts(:,:,ii);
poiss = lambda_p(:,:,ii);
identity = lambda_i(:,:,ii);


realGiniCoeff = ginicoeff(real, 1); %% across units
unitCycleGiniCoeff = ginicoeff(unitCycle, 1); 
timeSwapGiniCoeff = ginicoeff(timeSwap, 1);
poissGiniCoeff = ginicoeff(poiss, 1);
identityGiniCoeff = ginicoeff(identity, 1);



% histograms

temp = [realGiniCoeff unitCycleGiniCoeff timeSwapGiniCoeff poissGiniCoeff identityGiniCoeff];
bins = linspace(min(temp(:)), max(temp(:)), 50);

realH = hist(realGiniCoeff, bins);
UnitCycleH = hist(unitCycleGiniCoeff, bins);
timeSwapH = hist(timeSwapGiniCoeff, bins);
poissH = hist(poissGiniCoeff, bins);
identityH = hist(identityGiniCoeff, bins);


% statistical tests for comparing the distributions

realVScycleP = ranksum(realGiniCoeff, unitCycleGiniCoeff);
realVStimeSwapP = ranksum(realGiniCoeff, timeSwapGiniCoeff);
realVSpoissP = ranksum(realGiniCoeff, poissGiniCoeff);
realVSidentityP = ranksum(realGiniCoeff, identityGiniCoeff);


% smoothing and plot the distributions
sigma = 2;
halfwidth = 2 * sigma;
smoothwin = gausswindow(sigma, halfwidth);

realH = conv(realH, smoothwin, 'same');
UnitCycleH = conv(UnitCycleH, smoothwin, 'same');
timeSwapH = conv(timeSwapH, smoothwin, 'same');
poissH = conv(poissH, smoothwin, 'same');
identityH = conv(identityH, smoothwin, 'same');

% realH = realH ./ sum(realH);
% UnitCycleH = UnitCycleH ./ sum(UnitCycleH);
% timeSwapH = timeSwapH ./ sum(timeSwapH);
% poissH = poissH ./ sum(poissH);
% identityH = identityH ./ sum(identityH);

%%%%%%

hold on

plot(bins, realH, 'k', 'linewidth', 1.5)
plot(bins, UnitCycleH, 'c','linewidth', 1)
plot(bins, timeSwapH, 'b','linewidth', 1)
plot(bins, poissH, 'g','linewidth', 1)
plot(bins, identityH, 'm','linewidth', 1)


% significance bars
[realPy, realPx]= max(realH);
[UnitCyclePy, UnitCyclePx]= max(UnitCycleH);
[timeSwapPy, timeSwapPx]= max(timeSwapH);
[poissPy, poissPx]= max(poissH);
[identityPy, identityPx]= max(identityH);

firstmax = max([realPy, UnitCyclePy, timeSwapPy, poissPy, identityPy])*1.05;
textsize = 8;
% Wilcoxon ranksum test
if realVScycleP < 0.05
    y = firstmax;
    plot([bins(realPx) bins(UnitCyclePx)],[y y],'k','linewidth', 1.5)
    
    if realVScycleP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(UnitCyclePx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

end

if realVStimeSwapP < 0.05
    y = firstmax+0.5;
    plot([bins(realPx) bins(timeSwapPx)],[y y],'k','linewidth', 1.5)
    
    if realVStimeSwapP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(timeSwapPx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

end

if realVSpoissP < 0.05
    y = firstmax+1;
    plot([bins(realPx) bins(poissPx)],[y y],'k','linewidth', 1.5)
    
    if realVSpoissP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(poissPx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

    
end

if realVSidentityP < 0.05
    y = firstmax+1.5;
    plot([bins(realPx) bins(identityPx)],[y y],'k','linewidth', 1.5)
    
    if realVSidentityP < 0.01
        psign = '**';
    else
        psign = '*';
    end
        
    text(mean([bins(realPx), bins(identityPx)]), y+0.05, psign, 'fontsize', textsize, 'HorizontalAlignment', 'center')

end

hold off

set(gca, 'fontsize', 10)
xlabel('Gini coefficient', 'fontsize', 10)
ylabel('Number of states', 'fontsize', 10)
title('Lambda sparsity (across units)', 'fontsize', 10)
