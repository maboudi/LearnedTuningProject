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

% 
% % for RatU
% electrodePositions = [...
% -7.5    105    0; ...
% 7.5     97.5   0; ...
% -7.5    90     0; ...
% 7.5     82.5   0; ...
% -7.5    75     0; ...
% 7.5     67.5   0; ...
% -7.5    60     0; ...
% 7.5     52.5   0; ...
% -7.5    45     0; ...
% 7.5     37.5   0; ...
% -7.5    30     0; ...
% 7.5     22.5   0; ...
% -7.5    15     0; ...
% 7.5     7.5    0; ...
% -7.5    0      0; ...
% 7.5     -7.5   0; ...
% ];
% 
% % for one of the probes 8 shanks and for the other probe 4 shanks were
% % successfully recorded
% 
% nShanksPerProbe = 8;


p_base_shank = [0 -22; 14.5 10; 23.5 180; 23.5 300; -23.5 300; -23.5 180; -14.5 10; 0 -22];
p_base_elec  = [-5 -7.5; 5 -7.5; 5 7.5; -5 7.5; -5 -7.5];



% % for each unit we need the average spike waveform on each site

filename = 'Achilles_10252013.spikes.cellinfo.mat';
load(filename, 'spikes')

waveforms = spikes.indiv_waveforms;

nUnits = numel(spikes.UID);

nSpikes = size(waveforms{1}, 3);

minPnt = nan(nUnits, 3, nSpikes);

nElecs = size(electrodePositions, 1);

fprintf('\n estimate cell postion for each individual spike ')

for iu = 1: nUnits
    
    iu
    wavs = waveforms{iu}(:, :, 1:min(2000, size(waveforms{iu}, 3)));
    nSpikes = size(wavs, 3);
    
    temp = zeros(nSpikes, 3);
    parfor ispk = 1:nSpikes
        
%         if mod(ispk, 1000) == 0
%            fprintf('\b\b\b.%d%%', floor(ispk/nSpikes*100)) 
%         end
        
        currSpikeWaves = wavs(:, :, ispk);
    
        peaks = min(currSpikeWaves, [], 2);
        peaks = double(peaks);
    
        % gradient descent

        estPos = [0 0 -50]; % If we initiate z at zero it doesn't converge to the global minimum (compare the results with greedy search)
        alpha = 1e-5;
        nIters = 1000;
        [estPos, J_history] = gradientDecent(electrodePositions, peaks, estPos, alpha, nIters);
        
        temp(ispk, :) = estPos;
    end 
    
    minPnt(iu, :, 1:nSpikes) = permute(temp, [3 2 1]);
end


mean_minPnt = nanmean(minPnt, 3);
std_minPnt  = nanstd(minPnt, [], 3);


%%

nDots = 2000;

colors = {'g'; 'r'; 'b'; 'k'; 'm'};

figure;
hold on

for ihem = 1%:2
    
    
    ax(ihem) = subplot(1,1,ihem);
    
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
    
        
        % % %
%         spikePos = minPnt(spikes.shankID(1: nUnits) == currShankID, :, :);
%         spikePos(:, 1, :) = spikePos(:, 1, :)+(iShank-1)*200;
%         
%         for iUnit = 1:5
%             
%            currDots = permute(squeeze(spikePos(iUnit, :, :)), [2, 1]); 
%            scatter3(currDots(1:nDots, 1), currDots(1:nDots, 2), currDots(1:nDots, 3), 1, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor', colors{iUnit})
%            
%         end
        

        cellPos = mean_minPnt(spikes.shankID(1: nUnits) == currShankID, 1:2);
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


cl = colormap('jet');
maxZz = -min(mean_minPnt(:, 3));

figure;


for ii = 1:30
    
    UnitID = randi(nUnits, 1);
    
    subplot(3,10, ii)

    cellPos = mean_minPnt(UnitID, :);
    
    if cellPos(3) == 0; cellPos(3) = -1; end
    try
    scatter(cellPos(1), cellPos(2), 50, 'MarkerFaceColor', cl(floor(-cellPos(3)/maxZz * 64)+1, :), 'MarkerEdgeColor', 'none')
    catch
        continue
    end
    
%     wavs = waveforms{UnitID};
    
    wavs = mean(spikes.indiv_waveforms{UnitID}, 3);
    
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

    

    for jj = 1:nElecs

        elecP = electrodePositions(jj, :);

       if mod(jj, 2) == 1
           plot(elecP(1)-30-15: elecP(1)-30+16, wavs(jj, :)*20/maxApm+elecP(2), 'k', 'linewidth', 3)
       else
           plot(elecP(1)+30-15: elecP(1)+30+16, wavs(jj, :)*20/maxApm+elecP(2), 'k', 'linewidth', 3)
       end


    end
    
    xlim([-100 100])
    
    
    
%      axis equal
     axis off
     title(sprintf('unit %d', UnitID), 'fontsize', 12)
    
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



function [estPos, J_history] = gradientDecent(elecPositions, wavePeakAmps, estPos, alpha, nIters)


% alpha = 0.01;
% nIters  = 1500;
epsilon = 1;

nElecs = size(elecPositions, 1);

J_history = nan(nIters, 1);
diverged  = 0;
iIter = 1;

while iIter < nIters && ~diverged
    
    
    
    diffx = nan(nElecs^2/2 - nElecs, 1);
    diffy = nan(nElecs^2/2 - nElecs, 1);
    diffz = nan(nElecs^2/2 - nElecs, 1);
    
    count = 1;
   for ii = 1:nElecs
       for jj = (ii+1) : nElecs
           
           temp = wavePeakAmps(ii)*sum((estPos - elecPositions(ii,:)).^2) - wavePeakAmps(jj)*sum((estPos - elecPositions(jj,:)).^2);
           if temp >= 0
               diffx(count) = (wavePeakAmps(ii)*(estPos(1) - elecPositions(ii,1))) - (wavePeakAmps(jj)*(estPos(1) - elecPositions(jj,1)));
               diffy(count) = (wavePeakAmps(ii)*(estPos(2) - elecPositions(ii,2))) - (wavePeakAmps(jj)*(estPos(2) - elecPositions(jj,2)));
               diffz(count) = (wavePeakAmps(ii)*(estPos(3) - elecPositions(ii,3))) - (wavePeakAmps(jj)*(estPos(3) - elecPositions(jj,3)));
               
           elseif temp < 0
               diffx(count) = -(wavePeakAmps(ii)*(estPos(1) - elecPositions(ii,1))) + (wavePeakAmps(jj)*(estPos(1) - elecPositions(jj,1)));
               diffy(count) = -(wavePeakAmps(ii)*(estPos(2) - elecPositions(ii,2))) + (wavePeakAmps(jj)*(estPos(2) - elecPositions(jj,2)));
               diffz(count) = -(wavePeakAmps(ii)*(estPos(3) - elecPositions(ii,3))) + (wavePeakAmps(jj)*(estPos(3) - elecPositions(jj,3)));
           end
               
           count = count +1;
       
       end
   end
   
   
   estPos_old = estPos;
   
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


