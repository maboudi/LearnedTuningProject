% [PlaceMap, OccupancyMap] = PFClassic(Pos, SpkCnt, Smooth, nGrid, TopRate)
% Place field calculation "classic style" where a smoothed spike count map
% is divided by a smoothed occupancy map
%
% Pos is a nx2 array giving position in each epoch.  It should be in
% the range 0 to 1.  SpkCnt gives the number of spikes in each epoch.
%
% Smooth is the width of the Gaussian smoother to use (in 0...1 units).
%
% nGrid gives the grid spacing to evaluate on (should be larger than 1/Smooth)
%
% TopRate is for the maximum firing rate on the color map (if you display it)
% if you don't specify this it will be the maximum value
%
% optional output OccupancyMap is a smoothed occupancy map

function [PlaceMap, sTimeSpent] = PFClassic(Pos, SpkCnt, Smooth, nGrid, TopRate)
% 
% if (max(Pos(:,1)>1));
%     min1 = min(Pos(:,1)); max1 = max(Pos(:,1));
%     min2 = min(Pos(:,2)); max2 = max(Pos(:,2));
%     min1 = 0; min2 = 0;
%     Pos(:,1) = (Pos(:,1)-min1)/(max1-min1);
%     Pos(:,2) = (Pos(:,2)-min2)/(max2-min2);
% end
%     

% integrized Pos (in the range 1...nGrid
iPos = 1+floor(nGrid*Pos/(1+eps));

% make unsmoothed arrays
TimeSpent = full(sparse(iPos(:,1), iPos(:,2), 1, nGrid, nGrid));
nSpikes = full(sparse(iPos(:,1), iPos(:,2), SpkCnt, nGrid, nGrid));

% do the smoothing
r = (-nGrid:nGrid)/nGrid;

% Smooth = .001;
Smoother = exp(-r.^2/Smooth^2/2);
% 
sTimeSpent = conv2(Smoother, Smoother, TimeSpent, 'same');

% Smooth = .002;
% Smoother = exp(-r.^2/Smooth^2/2);
snSpikes = conv2(Smoother, Smoother, nSpikes, 'same');
% 
% sTimeSpent = full(sparse(iPos(:,1), iPos(:,2), 1, nGrid, nGrid));
% snSpikes = full(sparse(iPos(:,1), iPos(:,2), SpkCnt, nGrid, nGrid));

PlaceMap = snSpikes./(sTimeSpent+eps);
% NB not regularized to be mean firing rate in non-visited areas.

if nargout==0
    FireRate = PlaceMap*32552/512;
	if nargin<5 
        TopRate = [];
	end
	PFPlot(PlaceMap, sTimeSpent, TopRate, 1);
end

