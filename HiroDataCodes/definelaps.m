function [LRlaptimings, RLlaptimings] = definelaps(position, time, inboundPos, FileBase)



posInd = zeros(length(time), 1); %% assigning an index based on the location of the position points relative to the track ends

posInd(find(position < inboundPos(1))) = 1; %% if within trach end 1 then the index is 1 
posInd(find(position > inboundPos(2))) = 2; 
                                            %%% if neither of the ends (could be on the track or rest box assign 0)

 
                                            
passInd = [0; diff(posInd)];   % type of transition at the positions
                                 % end1 outbound --> -1
                                 % end1 inbound  --> +1
                                 % end2 outbound --> -2
                                 % end2 inbound  --> +2

passes = find(passInd); %% points with non-zero passInd



passTypes = passInd(passes); %% the type of pass at the corresponding positions


                                             
                                             
LRstart = find(passTypes(1:end-1) == -1 & passTypes(2:end) == +2); %% the start positions of Left(end1) to right(end2) traversals 
LRsegments = [passes(LRstart) passes(LRstart+1)]; %% the start and end points in each segment


RLstart = find(passTypes(1:end-1) == -2 & passTypes(2:end) == +1); %% the start positions of right(end2) to left(end1) traversals 
RLsegments = [passes(RLstart) passes(RLstart+1)]; %% the start and end points in each segment


%%% smoothing the position 

sigma = 5; %% do a little smoothing
halfwidth = 15;

win = gausswindow(sigma, halfwidth);
position = conv(position, win, 'same');
LRlaps = zeros(size(LRsegments));
RLlaps = zeros(size(RLsegments));

for lap =1 : size(LRsegments, 1)
    
    startpoint = LRsegments(lap, 1);
    endpoint = LRsegments(lap, 2);
    
    
    precedpoints = startpoint - 1000 : startpoint;
    possdiff = [0 diff(position(precedpoints))];
    maxInd = find(possdiff .* [possdiff(2:end) 0] < 0, 1, 'last');
    
    if ~isempty(maxInd)
        LRlaps(lap, 1) = precedpoints(maxInd);
    else
        LRlaps(lap, 1) = precedpoints(1);
    end
     
    
    postpoints = endpoint : endpoint + 1000;
    possdiff = [0 diff(position(postpoints))];
    maxInd = find(possdiff .* [possdiff(2:end) 0] < 0, 1, 'first');
    
    if ~isempty(maxInd)
        LRlaps(lap, 2) = postpoints(maxInd);
    else
        LRlaps(lap, 2) = postpoints(end);
    end
end


for lap = 1 : size(RLsegments, 1)
    
    startpoint = RLsegments(lap, 1);
    endpoint = RLsegments(lap, 2);
    
    
    precedpoints = startpoint - 1000 : startpoint;
    possdiff = [0 diff(position(precedpoints))];
    maxInd = find(possdiff .* [possdiff(2:end) 0] < 0, 1, 'last');
    
    if ~isempty(maxInd)
        RLlaps(lap, 1) = precedpoints(maxInd);
    else
        RLlaps(lap, 1) = precedpoints(1);
    end

    postpoints = endpoint : endpoint + 1000;
    possdiff = [0 diff(position(postpoints))];
    maxInd = find(possdiff .* [possdiff(2:end) 0] < 0, 1, 'first');
    
    if ~isempty(maxInd)
        RLlaps(lap, 2) = postpoints(maxInd);
    else
        RLlaps(lap, 2) = postpoints(end);
    end
    
end

% LRlaps = LRlaps(find(LRlaps(:,1) .* LRlaps(:,2)), :);
% RLlaps = RLlaps(find(RLlaps(:,1) .* RLlaps(:,2)), :);
 

figure;
hold on
plot(time/1e6, position, 'k', 'linewidth', 2)

for ii = 1:length(LRlaps)
    
    plot(time(LRlaps(ii, 1):LRlaps(ii, 2))/1e6, position(LRlaps(ii, 1):LRlaps(ii, 2)), 'b', 'linewidth', 3)
    
end

for ii = 1:length(RLlaps)
    plot(time(RLlaps(ii, 1):RLlaps(ii, 2))/1e6, position(RLlaps(ii, 1):RLlaps(ii, 2)), 'r', 'linewidth', 3)
end

hold off
set(gca, 'fontsize', 16)
xlabel('time(sec)', 'fontsize', 20)
ylabel('position', 'fontsize', 20)

saveas(gcf, [FileBase '/laps.fig'])


LRlaptimings = time(LRlaps);
RLlaptimings = time(RLlaps);

end









