IC_set1 = ICweights.coarseTS';
IC_set2 = ICweights.fineTS';

similarity = zeros(size(IC_set1, 2), size(IC_set2, 2));
for ii = 1:size(IC_set1, 2)
    for jj = 1:size(IC_set2, 2)
        
        similarity(ii, jj) = sum(IC_set1(:, ii) .* IC_set2(:, jj))/norm(IC_set1(:, ii))/norm(IC_set2(:, jj));
        
    end
end

cm = redblue;
figure; 

subplot(1,3,1)
imagesc(IC_set1); colormap('jet')
xlabel('ICs', 'fontsize', 12)
ylabel('units', 'fontsize', 12)
title('125 ms ICs', 'fontsize', 12)

subplot(1,3,2)
imagesc(IC_set2); colormap('jet')
xlabel('ICs', 'fontsize', 12)
ylabel('units', 'fontsize', 12)
title('20 ms ICs', 'fontsize', 12)

subplot(1,3,3)
imagesc(similarity); colormap(cm)
xlabel('125 ms ICs', 'fontsize', 12)
ylabel('20 ms ICs', 'fontsize', 12)
title('cosine similarity', 'fontsize', 12)




ICs_actual    = ICweights.coarseTS';
ICs_surrogate = ICweights_s.coarseTS';

ICs_actual    = ICs_actual .* sign(ICs_actual);
ICs_surrogate = ICs_surrogate .* sign(ICs_surrogate);

giniCoeff_actual    = ginicoeff(ICs_actual);
giniCoeff_surrogate = ginicoeff(ICs_surrogate);


median(giniCoeff_actual)
median(giniCoeff_surrogate)    