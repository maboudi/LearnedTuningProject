% clear;
% clc;

ii = 100;

load('actual.mat', 'lambdas', 'transmats')
Acttransmat = transmats(:,:,ii);
lambda = lambdas(:,:,ii);

ActualPnts = MultiDScale(lambda);
[ActLen, Actualseq] = pathLen_thresh2(Acttransmat, 0.07);
Actualseq = cell2mat(Actualseq);   

load('poisson.mat', 'lambdas', 'transmats')
Potransmat = transmats(:,:,ii);
lambda = lambdas(:,:,ii);

PoissonPnts = MultiDScale(lambda);
[PoissLen, Poissonseq] = pathLen_thresh2(Potransmat, 0.15);
Poissonseq = cell2mat(Poissonseq); 

load('timeCycle.mat', 'lambdas', 'transmats')
Cytransmat = transmats(:,:,ii);
lambda = lambdas(:,:,ii);

CyclePnts = MultiDScale(lambda);
[CycleLen, Cycleseq] = pathLen_thresh2(Cytransmat, 0.15);
Cycleseq = cell2mat(Cycleseq); 

load('timeswap.mat', 'lambdas', 'transmats')
Swtransmat = transmats(:,:,ii);
lambda = lambdas(:,:,ii);

SwapPnts = MultiDScale(lambda);
[SwapLen, Swapseq] = pathLen_thresh2(Swtransmat, 0.08);
Swapseq = cell2mat(Swapseq); 


%% Plot the longest path through the coactivity space 

colors = [[50 150 50]/255; [150 0 0]/255; [30 144 200]/255; [180 150 110]/255];


figure;
hold on
%actual
plot3(ActualPnts(:,1), ActualPnts(:,2), ActualPnts(:,3), 'color', colors(1,:), 'marker', 'o', 'LineStyle', 'none', 'markersize', 7)
h1 = plot3(ActualPnts(Actualseq, 1), ActualPnts(Actualseq, 2), ActualPnts(Actualseq, 3), 'color', colors(1,:), 'linewidth', 1);

%poisson
plot3(PoissonPnts(:,1), PoissonPnts(:,2), PoissonPnts(:,3),'color', colors(4,:), 'marker', 'square', 'LineStyle', 'none', 'markersize', 7)
h2 = plot3(PoissonPnts(Poissonseq, 1), PoissonPnts(Poissonseq, 2), PoissonPnts(Poissonseq, 3),'color', colors(4,:), 'linewidth', 1);

%temporal shuffle
plot3(CyclePnts(:,1), CyclePnts(:,2), CyclePnts(:,3),'color', colors(3,:), 'marker', 'square', 'LineStyle', 'none', 'markersize', 7)
h4 = plot3(CyclePnts(Cycleseq, 1), CyclePnts(Cycleseq, 2), CyclePnts(Cycleseq, 3),'color', colors(3,:), 'linewidth', 1);


%time swap
plot3(SwapPnts(:,1), SwapPnts(:,2), SwapPnts(:,3), 'color', colors(2,:), 'marker', '^', 'LineStyle', 'none', 'markersize', 7)
h3 = plot3(SwapPnts(Swapseq, 1), SwapPnts(Swapseq, 2), SwapPnts(Swapseq, 3), 'color', colors(2,:), 'linewidth', 1);

xlabel('dimension 1', 'fontsize', 16)
ylabel('dimention 2', 'fontsize', 16)
zlabel('dimention 3', 'fontsize', 16)

legend([h1, h2, h3, h4], {'actual', 'Poisson', 'time swap', 'temporal shuffle'});



%% Plot all the transitions 


% actual
figure;
hold on
plot3(ActualPnts(:,1), ActualPnts(:,2), ActualPnts(:,3), 'color', colors(1,:), 'marker', 'o', 'LineStyle', 'none', 'markersize', 7)


Acttransmat(Acttransmat<0.1) = 0;

noStates = size(Acttransmat);

for ii = 1:noStates
    for jj = 1:noStates
        
        if jj ~= ii & Acttransmat(ii, jj) > 0
            plot3(ActualPnts([ii jj], 1), ActualPnts([ii jj], 2), ActualPnts([ii jj], 3), 'linewidth', Acttransmat(ii, jj)*10, 'color', colors(1,:))
    
        end
    end
end




%% Functions
function points = MultiDScale(lambda)

[noUnits, noStates] = size(lambda);

similarity = zeros(noStates);
expansionParameter = noUnits;

for ii = 1 : noStates
    
    vector1 = lambda(:, ii);
    vector1 = vector1/sum(vector1);
    
    for jj = 1 : noStates
        
        vector2 = lambda(:, jj);
        vector2 = vector2/sum(vector2);
        
        if jj ~= ii
            similarity(ii, jj) = (1+vector1'*vector2)^expansionParameter;
        end
        
    end   
end
distance = 1./similarity;
distance(eye(noStates)== 1) = 0;

points = mdscale(distance, 3);
end




    