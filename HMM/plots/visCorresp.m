function visCorresp(binnedfRates, M1M2corresp, pbeStateProbM1, pbeStateProbM2, noPBE2plot) 


nPBEs = size(binnedfRates, 1); 

rndselectIdx = randi(nPBEs, noPBE2plot, 1);
rndselectPBEs = binnedfRates(rndselectIdx, 1); % one milisecond bins for raster plot

rndselectStateProbM1 = pbeStateProbM1(rndselectIdx);
rndselectStateProbM2 = pbeStateProbM2(rndselectIdx);


pbeLen = zeros(noPBE2plot, 1);
for pbe = 1: noPBE2plot
    pbeLen(pbe) = size(binnedfRates{rndselectIdx(pbe), 2}, 2);
end 


rndselectPBEs = cell2mat(rndselectPBEs');

rndselectStateProbM1 = cell2mat(rndselectStateProbM1');
rndselectStateProbM2 = cell2mat(rndselectStateProbM2');

nStatesM1 = size(rndselectStateProbM1, 1);
nStatesM2 = size(rndselectStateProbM2, 1);

rndselectStateProbM1 = rndselectStateProbM1 ./ repmat(sum(rndselectStateProbM1, 1), [nStatesM1 1]); % in case we are using observation probabilities instead of gamma distibutions as state probability distribution
rndselectStateProbM2 = rndselectStateProbM2 ./ repmat(sum(rndselectStateProbM2, 1), [nStatesM2 1]);

[~, binMostProbM1State] = max(rndselectStateProbM1);
[~, binMostProbM2State] = max(rndselectStateProbM2);


nUnits = size(rndselectPBEs, 1);
nSelectTBins = size(rndselectStateProbM1, 2);

sigBin = zeros(nSelectTBins, 1);
binSlefTransM1 = zeros(nSelectTBins, 1);
binSlefTransM2 = zeros(nSelectTBins, 1);

for bin = 1: nSelectTBins
   if M1M2corresp(binMostProbM1State(bin), binMostProbM2State(bin)) > -log10(0.05/nStatesM2/nStatesM1)
      
       sigBin(bin) = 1;
   end
   
    if bin > 1
        if binMostProbM1State(bin) == binMostProbM1State(bin-1)
            binSlefTransM1(bin) = 1;
        end
        
        if binMostProbM2State(bin) == binMostProbM2State(bin-1)
            binSlefTransM2(bin) = 1;
        end 
    end
end


cumpbeLen = cumsum(pbeLen);


figure; 

%
h(1) = subplot(2,1,1);  
rndselectPBEs2 = rndselectPBEs + repmat((1:nUnits)', [1, size(rndselectPBEs, 2)]);
rndselectPBEs2(rndselectPBEs == 0) = 0;

imagesc(rndselectPBEs2)
hold on
plot([(1:nSelectTBins-1)*20+0.5; (1:nSelectTBins-1)*20+0.5], repmat([0 nUnits]', 1, nSelectTBins-1), '--', 'color', [0.7 0.7 0.7], 'linewidth', 1)
hold on
plot([cumpbeLen(1:end-1)*20+0.5 cumpbeLen(1:end-1)*20+0.5]', repmat([0 nUnits]', 1, noPBE2plot-1), '-w', 'linewidth', 2)

mymap = colormap('jet');
mymap = mymap(randperm(length(colormap)), :);
mymap(1, :) = [0.1 0.1 0.1];

colormap(h(1), mymap)

%
h(2) = subplot(2,1,2);
imagesc((0:nSelectTBins-1)*20+10+0.5, 1:(nStatesM2+nStatesM1),[rndselectStateProbM1; rndselectStateProbM2])
hold on
plot([(1:nSelectTBins-1)*20+0.5; (1:nSelectTBins-1)*20+0.5], repmat([0 nStatesM2+nStatesM1]', 1, nSelectTBins-1), '--', 'color', [0.7 0.7 0.7], 'linewidth', 1)
hold on
plot([cumpbeLen(1:end-1)*20+0.5 cumpbeLen(1:end-1)*20+0.5]', repmat([0 nStatesM2+nStatesM1]', 1, noPBE2plot-1), '-k', 'linewidth', 2); 

c = gray;
c = flipud(c);
colormap(h(2), c);

line([0 nSelectTBins*20], [nStatesM1-0.5 nStatesM1-0.5], 'color', 'k', 'linewidth', 2)

for bin = 1: nSelectTBins
   if sigBin(bin)
       p = patch([(bin-1)*20+0.5 bin*20+0.5 bin*20+0.5 (bin-1)*20+0.5],[0 0 nStatesM2+nStatesM1 nStatesM2+nStatesM1], ...
            'r','edgecolor', 'none');
       set(p,'FaceAlpha', 0.1)
   end
   
   if binSlefTransM1(bin)
      line([(bin-1)*20+0.5 bin*20+0.5], [nStatesM1-5 nStatesM1-5], 'color', 'g', 'linewidth', 3)
   end
   
   if binSlefTransM2(bin)
       line([(bin-1)*20+0.5 bin*20+0.5], [nStatesM1+nStatesM2-5 nStatesM1+nStatesM2-5], 'color', 'b', 'linewidth', 3)
   end
  
end

linkaxes(h, 'x')

end
