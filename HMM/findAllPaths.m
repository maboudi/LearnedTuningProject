function [visitedNodes, paths] = findAllPaths(startNodes, transmat, visitedNodes, paths)

for i = 1 : length(startNodes)
    
    visitedNodes = [visitedNodes startNodes(i)];
    connections = find(transmat(visitedNodes(end), :)); 
    connections = setdiff(connections, visitedNodes); %% excluding the already visited nodes
    
    if ~isempty(connections)
        [visitedNodes, paths] = findAllPaths(connections, transmat, visitedNodes, paths);
    else
        
        paths = [paths; visitedNodes];
    end
    visitedNodes(end) = [];
end

end