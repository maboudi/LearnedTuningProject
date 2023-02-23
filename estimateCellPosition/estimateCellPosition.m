electrodePositions = [...
-18.5  140  0; ...
16.5   120  0; ...
-14.5  100  0; ...
12.5   80   0; ...
-10.5  60   0; ...
8.5    40   0; ...
-8.5   20   0; ...
0      0    0]; % recording site positons for probe Buzsaki64

nShanksPerProbe = 8;


electrodePositions = [...
-16    180    0; ... 
14.75  160    0; ...
-13.5  140    0; ...
12.25  120    0; ...
-11    100    0; ...
9.75   80     0; ...
-8.5   60     0; ...
7.25   40     0; ...
-7    20      0; ...
0      0      0];

nShanksPerProbe = 6;


p_base_shank = [0 -22; 14.5 10; 23.5 180; 23.5 300; -23.5 300; -23.5 180; -14.5 10; 0 -22];
p_base_elec  = [-5 -7.5; 5 -7.5; 5 7.5; -5 7.5; -5 -7.5];



% % for each unit we need the average spike waveform on each site

filename = 'Achilles_10252013.spikes.cellinfo.mat';
load(filename, 'spikes')

waveforms = spikes.rawWaveform2;

nUnits = numel(spikes.UID);

% for greedy search
xrange = -200:2:200;
yrange = -200:2:(electrodePositions(1, 2)+ 200);
zrange = -200:2:200;


minPnt = nan(nUnits, 3);
nElecs = size(electrodePositions, 1);
tic
for iu = 1: nUnits
    
    wavs = waveforms{iu};
    
    peaks = min(wavs, [], 2);

    
    % gradient descent
    
    estPos = [0 0 -50]; % If we initiate z at zero it doesn't converge to the global minimum (compare the results with greedy search)
    alpha = 1e-5;
    nIters = 1000;
    [estPos, J_history] = gradientDecent(electrodePositions, peaks, estPos, alpha, nIters);
    minPnt(iu, :) = estPos;
    
    
    
    % greedy search
%     
%     sumDiff = nan(numel(xrange), numel(yrange), numel(zrange));
%     for ix = 1:numel(xrange)
%         
%         currX = xrange(ix);
%         
%         for iy = 1:numel(yrange)
%             currY = yrange(iy);
%             
%             for iz = 1:numel(zrange)
%                 currZ = zrange(iz);
%                 
%                 distSq = nan(nElecs, 1);
%                 for ielec = 1: size(electrodePositions, 1)
%                     distSq(ielec) = (currX - electrodePositions(ielec,1))^2 + (currY - electrodePositions(ielec,2))^2 + (currZ - electrodePositions(ielec,3))^2;
%                 end
%                 
%                 ampDistMulti = distSq .* peaks;
%                 
%                 diffs = zeros(nElecs^2/2 - nElecs, 1);
%                 count = 1;
%                 for ii = 1:nElecs
%                    
%                     for jj = (ii+1) : nElecs
%                        
%                        diffs(count) = abs(ampDistMulti(ii) - ampDistMulti(jj));
%                        count = count + 1;
%                     end
%                 end
%                 sumDiff(ix, iy, iz) = sum(diffs);
%                 
%             end
%         end
%     end
%     
%     [v,loc] = min(sumDiff(:));
%     [ixMax,iyMax, izMax] = ind2sub(size(sumDiff),loc);
%     
%     minPnt(iu, :) = [xrange(ixMax) yrange(iyMax) zrange(izMax)];  
%     
end
toc

%%
figure;
hold on

for ihem = 1:2
    
    ax(ihem) = subplot(2,1,ihem);
    hold on
    for iShank = 1:nShanksPerProbe

        p = p_base_shank;
        p(:, 1) = p(:, 1) + (iShank-1)*200;
        pgon = polyshape(p);

        plot(pgon, 'FaceColor', [0.3 0.3 .3], 'EdgeColor', 'none')

        for ielec = 1:nElecs

            p = p_base_elec + repmat(electrodePositions(ielec, 1:2), [5, 1]);
            p(:, 1) = p(:, 1) + (iShank-1)*200;
            pgon = polyshape(p);

            plot(pgon, 'FaceColor', [1 1 1], 'EdgeColor', 'k')


        end
        
        if nShanksPerProbe == 6
            currShankID = (ihem-1)*7 + iShank;
        else
            currShankID = (ihem-1)*8 + iShank;
        end
        
        cellPos = minPnt(spikes.shankID == currShankID, 1:2);
        cellPos(:,1) = cellPos(:,1)+(iShank-1)*200;


        scatter(cellPos(:, 1), cellPos(:, 2), 20 , 'filled')


    end
    axis equal
    axis off
end   

linkaxes(ax, 'xy')

%%
% % for an example unit draw waveforms next to corresponding electrodes on
% the shank

plotwidth  = 800;
plotheight = 500;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



cl = colormap('jet');
maxZz = -min(minPnt(:, 3));

for ii = 1:30
    
    UnitID = randi(nUnits, 1);
    
    subplot(3,10, ii)

    cellPos = minPnt(UnitID, :);
    wavs = waveforms{UnitID};
    maxApm = -min(wavs(:));
    
    
    hold on

    p = p_base_shank;
    pgon = polyshape(p);

    plot(pgon, 'FaceColor', [0.3 0.3 .3], 'EdgeColor', 'none')


    for ielec = 1:nElecs

        p = p_base_elec + repmat(electrodePositions(ielec, 1:2), [5, 1]);
        pgon = polyshape(p);

        plot(pgon, 'FaceColor', [1 1 1], 'EdgeColor', 'k')
    end

    if cellPos(3) == 0; cellPos(3) = -1; end
    scatter(cellPos(1), cellPos(2), -cellPos(3)+3, 'MarkerFaceColor', cl(floor(-cellPos(3)/maxZz * 64)+1, :), 'MarkerEdgeColor', 'none')


    for jj = 1:nElecs

        elecP = electrodePositions(jj, :);

       if mod(jj, 2) == 1
           plot(elecP(1)-30-15: elecP(1)-30+16, wavs(jj, :)*20/maxApm+elecP(2), 'k', 'linewidth', 1)
       else
           plot(elecP(1)+30-15: elecP(1)+30+16, wavs(jj, :)*20/maxApm+elecP(2), 'k', 'linewidth', 1)
       end


    end
    
    xlim([-100 100])
    
%      axis equal
     axis off
     title(sprintf('unit %d', UnitID), 'fontsize', fontsize)
    
end

%%

figure;
hold on

imagesc(xrange, yrange, mean(sumDiff, 3))
colormap('jet')


p = p_base_shank;
pgon = polyshape(p);

plot(pgon, 'FaceColor', [0.3 0.3 .3], 'EdgeColor', 'none')


for ielec = 1:nElecs

    p = p_base_elec + repmat(electrodePositions(ielec, 1:2), [5, 1]);
    pgon = polyshape(p);

    plot(pgon, 'FaceColor', [1 1 1], 'EdgeColor', 'k')
end
axis off
axis equal

UnitID = 4;
cellPos = minPnt(UnitID, :);
wavs = waveforms{UnitID};
maxApm = -min(wavs(:));


for jj = 1:nElecs

    elecP = electrodePositions(jj, :);

   if mod(jj, 2) == 1
       plot(elecP(1)-30-15: elecP(1)-30+16, wavs(jj, :)*20/maxApm+elecP(2), 'k', 'linewidth', 3)
   else
       plot(elecP(1)+30-15: elecP(1)+30+16, wavs(jj, :)*20/maxApm+elecP(2), 'k', 'linewidth', 3)
   end


end



%% functions

function [estPos, J_history] = gradientDecent(elecPositions, wavePeakAmps, estPos, alpha, nIters)

% This function estimates cell position based on the spike amplitudes
% recorded on different channels of a silicon probe


epsilon = 1;

nElecs = size(elecPositions, 1);
nAllElecPairs = (nElecs^2 - nElecs)/2;


J_history = nan(nIters, 1);
diverged  = 0;
iIter = 1;

while iIter < nIters && ~diverged

    diffx = nan(nAllElecPairs, 1);
    diffy = nan(nAllElecPairs, 1);
    diffz = nan(nAllElecPairs, 1);

    diff_d = nan(nAllElecPairs, 3);
    
    count = 1;
   for ii = 1:nElecs
    
       d_ii = estPos - elecPositions(ii,:);

       for jj = (ii+1) : nElecs

           d_jj = estPos - elecPositions(jj,:);
           
           costFunction = wavePeakAmps(ii)*sum(d_ii.^2) - wavePeakAmps(jj)*sum(d_jj.^2);

           diff_d(count, :) = sign(costFunction)*(wavePeakAmps(ii)*d_ii - wavePeakAmps(jj)*d_jj);
               
           count = count+1;
       
       end
   end
   
   estPos = estPos - alpha*nansum(diff_d)
   
   estPos(1) = estPos(1) - alpha * nansum(diffx);
   estPos(2) = estPos(2) - alpha * nansum(diffy);
   estPos(3) = estPos(3) - alpha * nansum(diffz);
    
   J_history(iIter) = computeCost(estPos, elecPositions, wavePeakAmps);
   
   if iIter > 1
      if  abs(J_history(iIter)-J_history(iIter-1)) < epsilon
          diverged = 1;
      end
   end
   
   iIter = iIter + 1;
   
end

end


function value = computeCost(estPos, elecPositions, wavePeakAmps)

nElecs = size(elecPositions, 1);


distSq = nan(nElecs, 1);
for ielec = 1: nElecs
    distSq(ielec) = (estPos(1) - elecPositions(ielec,1))^2 + (estPos(2) - elecPositions(ielec,2))^2 + (estPos(3) - elecPositions(ielec,3))^2;
end

ampDistMulti = distSq .* wavePeakAmps;



diffs = zeros(nElecs^2/2 - nElecs, 1);
count = 1;

for ii = 1:nElecs
    for jj = (ii+1) : nElecs

       diffs(count) = abs(ampDistMulti(ii) - ampDistMulti(jj));
       count = count + 1;
    end
end

value = sum(diffs);


end


