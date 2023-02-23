function [unique_paths, pathLen] = pathFinder(transmat, transProbThresh, overlapThresh)



% developing an adjacency matrix for corresponding graph (pairs are connected if the the transtion probability between them is above a certain threshold)

transmat(transmat >= transProbThresh) = 1;
transmat(transmat <= transProbThresh) = 0;


% finding all the pathes in the adjacency matrix: staring from each node
% find all the sequences of nodes connected together. There is no cycle in
% any of the found paths

startNodes = 1:size(transmat, 1);
visitedNodes = [];
paths = cell(0);

[~, paths] = findAlPaths(startNodes, transmat, visitedNodes, paths);


% find the length of each path (number of nodes)

pathLen = zeros(length(paths), 1);
for ii = 1 : length(paths)
    pathLen(ii) = length(paths{ii});
end


rmvidx2 = find(pathLen < 2); % remove the paths(!) with just one nodes

pathLen(rmvidx2) = [];
paths(rmvidx2) = [];

unique_paths = paths;

% %%%%%%
% 
% 
% %%% find the paths which are parts of some others and eliminate them using
% %%% analysis of sequence
% 
% 
% % assigning each state a unique string, so different combinations of them
% % for the paths
% 
% alphabet = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j','x'};
% numbers = {'1', '2', '3', '4', '5', '6', '7', '8', '9'};
% 
% paths_strings = cell(length(paths), 1);
% for ip = 1 : length(paths)
%     temp = mat2str(paths{ip});
%     
%     temp = strrep(temp, '[', ' ');
%     temp = strrep(temp, ']', ' ');
%     
%     for i = 1:floor(size(transmat, 1)/10)
%         for j = 0:9
%             temp = strrep(temp, [num2str(i) num2str(j)], [alphabet{i} num2str(j)]);
%         end
%     end
%     
%     for ii = 1 : 9
%         temp = strrep(temp, [' ' num2str(ii) ' '], [' x' num2str(ii) ' ']);
%     end
%     
%     temp = strrep(temp, ' ', '');
%     paths_strings{ip} = temp;
%     
% end
% 
% 
% % for each path find a longer path contating it
% 
% isRedundant = zeros(length(paths), 1);
% overlappingPath = zeros(length(paths), 1);
% 
% rmvflag = zeros(1, length(paths));
% 
% for ip = 1:length(paths)
% %      ip
%     set2comp = setdiff(1:length(paths), [ip find(rmvflag)]); % the set of paths to compare: all the paths other than those already known as redundant (previous loops)
% 
%     str1 = paths_strings{ip};
%  
%     for ip2 = set2comp %% loop through the comparsion set
%         
%         str2 = paths_strings{ip2};
%         
%         if length(str2) >= length(str1) %% just compare with paths either of the same length or larger
%             substr = commonsubstring(str1, str2); %% find the longest common substring
% 
% 
%             if length(substr) > 1
%                 if ismember(substr(1), numbers)
%                    substr(1) = [];
%                 end
% 
%                 if ismember(substr(end), alphabet)
%                    substr(end) = [];
%                 end
% 
%                 overlap = length(substr)/length(str1); %% calulationg the ratio of sting1 overlapping with a larger path
% 
%             else
%                 overlap = 0;
%             end
% 
%             
%             % compare the measured overlap with a threshold
%             if overlap >= overlapThresh 
% 
% %                 overlap
% 
%                 isRedundant(ip) = overlap;
%                 overlappingPath(ip) = ip2; % the larger overlapping path
% 
%                 rmvflag(ip) = 1;
%                 
%                 break
%             end
%         else
%             continue
%         end
%         
%     end
%     
% end
% 
% % find the unique (non-redundant paths)
% unique_paths = paths;
% unique_paths(find(rmvflag)) = [];
% pathLen(find(rmvflag)) = [];
% 

end

