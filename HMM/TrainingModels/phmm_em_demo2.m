clear; clc;
O = 4; %% number of neurons in individual obseravtions
Data_bout = 11; %% number of time bins
Nobout = 200;
T = Nobout * Data_bout; 
Q = 6;  

lambda_curve = [...
    3 5 2 2 1 1 0 1 0 1 1; ...
    0 1 1 3 8 3 2 0 0 1 0; ...
    1 1 0 0 1 0 3 6 2 0 0; ...
    0 0 1 1 1 0 2 1 3 4 7];

lambda_state = [...
    4 2 0 1;...
    1 6 1 0;...
    1 0 5 0;...
    0 1 2 6];    

% lambda_curve = [4 0 0 0;0 2 0 0;0 0 3 0; 0 0 0 5];
% data_short = poissrnd(lambda_curve);
% data = repmat(data_short, [1 Nobout]);  
data = []; 
for ii = 1 : Nobout
     data = [data poissrnd(lambda_curve)];
 end
% initial guess of parameters
prior0 = normalise(rand(Q,1));
transmat0 = mk_stochastic(rand(Q,Q));

% Initialize lambda to a random data point
indices = randperm(T);

% lambda0 = data(:,indices(1:4))+ 0.1;
% lambda0 = [3 4 2 2; 0 3 1 5; 6 1 3 2;1 3 3 1];

lambda_state_1D = lambda_state(:);
for i = 1 : length(lambda_state_1D)
    
    if lambda_state_1D(i) ~= 0
        permutegen = randperm(ceil(lambda_state_1D(i)/2));
        rndgen = randn(1)/10;
    else
        permutegen = 1;
        rndgen = rand(1);
    end
    lambda0_1D(i) = lambda_state_1D(i) + sign(rndgen)*ceil(abs(rndgen))* permutegen(1);
    if lambda0_1D(i) == 0
        lambda0_1D(i) = 0.1;
    end
end
% lambda0 = reshape(lambda0_1D, size(lambda_state));

lambda0 = rand(4,6)*10;
 B = poisson_prob(data, lambda0);


[LL, prior1, transmat1, lambda1, nnn] = phmm_em(data, prior0, transmat0, lambda0, 'max_iter', 30);
size(find(floor(nnn)))

 B = poisson_prob(data, lambda1,1);
 path = viterbi_path(prior1, transmat1, B);
 
 %%% sorting the path and transmat
 
 unique_path(1) = path(1);
 j = 1;
 for i = 2 : length(path)
     if ~ismember(path(i), unique_path)
        j = j+1;
        unique_path(j) = path(i);
     end
 end
 sorted_unique_path = sort(unique_path);
 
 sorted_path = path;
 
 for j = 1 : Q
     sorted_path(find(path == unique_path(j))) = sorted_unique_path(j);
 end
 
 sorted_transmat = zeros(size(transmat1));
 
 for ii = 1 : Q
     for jj = 1 : Q
         sorted_transmat(ii, jj) = transmat1(unique_path(ii), unique_path(jj));
     end
 end
 
 sorted_lambda = zeros(size(lambda1));
 for kk = 1 : Q
     sorted_lambda(:,kk) = lambda1(:,unique_path(kk));
 end
 
 figure;
 subplot(4,2,[3 5])
 t = 1: length(sorted_path);
 plot(t(1:50),sorted_path(1:50),'--sk','linewidth',2)
 title('Most probable sequence of states','fontsize', 24)
 
 subplot(4,2,[4 6])
 colormap bone
 imagesc(sorted_transmat)
 title('States transition matrix','fontsize',24)

% loglik = mhmm_logprob(data, prior1, transmat1, mu1, Sigma1, mixmat1);
