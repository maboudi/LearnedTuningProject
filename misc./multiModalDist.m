function xThresh = multiModalDist(x, numComponents)

binWidth = 0.1;
binCtrs = min(x):binWidth:max(x);


counts = hist(x, binCtrs);
counts = counts./sum(counts);


if isempty(numComponents)

    AIC = zeros(1,5);
    obj = cell(1,5);

    for kk = 1:5
        obj{kk} = gmdistribution.fit(x, kk);
        AIC(kk) = obj{kk}.AIC;
    end

    [~, numComponents] = min(AIC);

end
numComponents

gaussMixParams = gmdistribution.fit(x, numComponents);


mu = [];
sigma = [];
pp = [];
for ii = 1:numComponents
    
    mu = [mu; gaussMixParams.mu(ii)];
    sigma = [sigma; gaussMixParams.Sigma(ii)];
    
    pp = [pp, gaussMixParams.PComponents(ii)];
    
end

[mu, sortIdx] = sort(mu, 'ascend');
sigma = sigma(sortIdx);
pp = pp(sortIdx);

sumgauss = gmdistribution(mu, permute(sigma, [3,2,1]), pp);


gauss1 = gmdistribution(mu(1), sigma(1));

%%%
% gauss2 = gmdistribution(mu(2), sigma(2));
% gauss3 = gmdistribution(mu(3), sigma(3));



figure;

sumgaussPdf = pdf(sumgauss, binCtrs');
sumgaussPdf = sumgaussPdf./sum(sumgaussPdf);

plot(binCtrs', sumgaussPdf, 'linewidth', 2, 'color', [.7 .7 .7])
hold on



lowxComp = pdf(gauss1, binCtrs');
lowxComp = lowxComp * max(counts)/max(lowxComp);
plot(binCtrs', lowxComp, 'r-', 'linewidth', 2)

hold on

% secondxComp = pdf(gauss2, binCtrs');
% secondxComp = secondxComp * max(counts)/max(secondxComp);
% plot(binCtrs', secondxComp, 'g-', 'linewidth', 2)
% 
% hold on
% 
% thirdxComp = pdf(gauss3, binCtrs');
% thirdxComp = thirdxComp * max(counts)/max(thirdxComp);
% plot(binCtrs', thirdxComp, 'b-', 'linewidth', 2)
% 
% hold on

plot(binCtrs, counts, 'o', 'markersize', 5, 'color', 'k')


% calculate the x threshold

xThresh = find(cumsum(lowxComp) < (1-1e-5) * sum(lowxComp), 1, 'last') * binWidth;


end
