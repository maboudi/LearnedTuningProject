function [mixingTime, bottleneckRatio, meanbRatio] = connectivity(transmat)

    
P = transmat;

numStates = size(P,1);

% Find the stationary distribution, which is the eigen vector corresponding
% to largest eigenvalue  of the transition matrix

[v,d] = eig(P');
eigenValues = diag(d); 
statDist = v(:, find(real(eigenValues) > 0.99))./sum(v(:, find(real(eigenValues) > 0.99)));

% Mixing time: How long does it take (or on average how many transitions are needed) to approximate enough the stationary
% distribution


maxtime = 40; 
dist = zeros(maxtime, 1);

for t = 1 : maxtime

    Pt = P^t; 
    dist(t) = 0.5*max(sum(abs(Pt - repmat(statDist', [numStates, 1])), 2));

end

mixingTime = find(dist < 0.25, 1, 'first');



% bottleneck ratio


pi = statDist;
sample = 0;
for iter = 1 : 50000
    
    S_size = randi(numStates-1);
    randgen = randperm(numStates);
    
    S = randgen(1:S_size);
    S_com = setdiff(1:numStates, S);
    
    Pi = sum(pi(S));
    
    if Pi > 0.5
        continue
    end
    
    sample = sample + 1;
    
    sumQ = 0;
    for i = 1:S_size
        for j = 1: length(S_com)
            sumQ = sumQ + pi(S(i)) * P(S(i), S_com(j));
        end
    end
    
    Phi(sample) = sumQ / Pi;
end

bottleneckRatio = min(Phi);
meanbRatio = median(Phi);
% bottleneckRatio = 1/(4 * mixingTime);

end


