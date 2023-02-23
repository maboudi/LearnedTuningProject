
sigma = 3; %% in milisecond as we had use 1 milisecond time bins
halfwidth = 5 * sigma;

smoothwin = gausswindow(sigma, halfwidth);


figure;

x0=0;
y0=0;
width=300;
height=300;
set(gcf,'units','points','position',[x0,y0,width,height])

hold on

dims = 10:10:80;
colors = colormap('jet');
counts_all = [];
group = [];

spatialInformation = cell(length(dims), 1);
for ii = 1 : 8
    
    load(['stateSpaceField_numStates = ' num2str(dims(ii))])
%     [counts, bins] = hist(spaceEntropy, 10);
%     
%     counts = conv(counts, smoothwin, 'same');
%     
%     counts = counts./sum(counts);
%     
%     plot(bins, counts, 'linewidth', 2, 'color', colors(ii*8,:))
    spatialInfo = zeros(1,dims(ii));
    for jj = 1:dims(ii)
        spatialInfo(jj) = - sum(posbin_Occ' .* gamma_sum2(jj, :) .* log2(gamma_sum2(jj, :)));
    end
    
    spatialInformation{ii} = spatialInfo;
    
%     counts = spatialInfo;
    counts = spaceEntropy';
    
    counts_all = [counts_all counts];
    group = [group ii*ones(1,length(counts))];
    
end






    