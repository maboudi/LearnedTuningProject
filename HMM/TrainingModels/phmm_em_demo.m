clear; clc;
O = 4; %% number of neurons in individual obseravtions
Data_bout = 11; %% number of time bins
Nobout = 20;
T = Nobout * Data_bout; 

%%%% initial test on a problem with few neurons and few states
Q = 23;  

% lambda_curve = [...
%     3 5 2 2 1 1 0 1 0 1 1; ...
%     0 1 1 3 8 3 2 0 0 1 0; ...
%     1 1 0 0 1 0 3 6 2 0 0; ...
%     0 0 1 1 1 0 2 1 3 4 7];

%%%%% simulation of 100 neurons

Noseg = 300;
neurons = 20;
init_lambda_curve = zeros(neurons, Noseg);

win_length = 100;
sigma = linspace(3,15,neurons);
sigma = sigma(randperm(neurons));

amp = linspace(5,10,neurons);
amp = amp(randperm(neurons));

center = [1:neurons] * round(Noseg/neurons);

for n = 1 : neurons
    win = amp(n) * gausswin(win_length, sigma(n));
    init_lambda_curve(n,center(n)) = 1;
    lambda_curve(n,:) = conv(win', init_lambda_curve(n,:));
end



 
data = []; 
for ii = 1 : Nobout
     data = [data poissrnd(lambda_curve)];
 end
% initial guess of parameters
prior0 = normalise(rand(Q,1));
transmat0 = mk_stochastic(rand(Q,Q));

% Initialize lambda to a random data point

lambda0 = rand(neurons,Q);
B = poisson_prob(data, lambda0);


[LL, prior1, transmat1, lambda1, nnn] = phmm_em(data, prior0, transmat0, lambda0, 'max_iter', 50);
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
 
 sorted_transmat2 = sorted_transmat(2:end, 2:end);
 for ii = 2 : Q-1
     sorted_transmat2(ii-1,ii) = sorted_transmat(ii,ii+1) + sorted_transmat(ii,1);
 end
 
 sorted_transmat2(Q-1,2) = sorted_transmat(Q,2) + sorted_transmat(Q,1);
 
 
 sorted_lambda = zeros(size(lambda1));
 for kk = 1 : Q
     sorted_lambda(:,kk) = lambda1(:,unique_path(kk));
 end
 
 %%%%%%
 
 indrnd = randperm(neurons);
%  figure;
%  for i = 1 : neurons
%     plot(1:size(sorted_lambda,2), sorted_lambda(indrnd(i),:),color(i))
%     hold on
%  end
%%%%%%
 
 figure;
 subplot(4,4,[1:4])

colormap jet
neuroncolor = colormap;
colorstep = floor(length(neuroncolor)/neurons);
for i = 1 : neurons
    plot(1:size(lambda_curve,2), lambda_curve(indrnd(i),:),'color', neuroncolor(i*colorstep,:))
    hold on
end
title('Neurons probability of firing as a function of time','fontsize', 20)
ylabel('Lambda','fontsize', 16)
% xlabel('Time','fontsize', 16)
 
 subplot(4,4,[5:8])
 for i = 1 : neurons
    plot(1:size(lambda_curve,2), lambda_curve(indrnd(i),:)*4,'color', neuroncolor(i*colorstep,:))
    hold on
 end
 t = 1: length(sorted_path);
 plot(t(1:400),sorted_path(1:400),'-','color', [0.6 0.6 0.6],'linewidth',1.5)
 title('Most probable sequence of states','fontsize', 16)
%  xlabel('Time','fontsize', 16)
 grid on
%  
%  subplot(4,4,[9 10 13 14])
%  colormap bone
%  imagesc(sorted_transmat)
%  title('States transition matrix','fontsize',16)
%  
%  subplot(4,4,[11 12 15 16])
%  imagesc(sorted_transmat2)
%  title('States transition matrix (modified)','fontsize',16)
 
  subplot(4,4,[10 11 14 15])
  figure;
  colormap bone
%   clim = [0 0.3]
 imagesc(sorted_transmat .* (ones(size(sorted_transmat,1)) - eye(size(sorted_transmat,1))))
 title('States transition matrix','fontsize',16)

% loglik = mhmm_logprob(data, prior1, transmat1, mu1, Sigma1, mixmat1);
