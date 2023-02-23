function [loglik, prior, transmat, lambda] = phmm_em(data, prior, transmat, lambda_old, varargin)

[max_iter, thresh, verbose, adj_prior, adj_trans, adj_lambda] = ...
    process_options(varargin, 'max_iter', 10, 'thresh', 1e-6, 'verbose', 1, ...
		    'adj_prior', 1, 'adj_trans', 1, 'adj_lambda', 1);
        
previous_loglik = -inf;
loglik = 0;
converged = 0;
num_iter = 1;
% LL = [];

if ~iscell(data)
  data = num2cell(data, [1 2]); % each elt of the 3rd dim gets its own cell
end
numex = length(data);


O = size(data{1},1);
Q = length(prior);

while (num_iter <= max_iter) & ~converged
  % E step
  
  [loglik, exp_num_trans, exp_num_visits1, lambda] = ...
      ess_phmm(prior, transmat, lambda_old, data);
  % M step
  
   if adj_prior
    prior = normalise(exp_num_visits1);
  end
  if adj_trans 
    transmat = mk_stochastic(exp_num_trans);
  end

  if adj_lambda
    
    lambda(find(floor(lambda*1000) == 0)) = 0.0001;
    lambda_old = lambda;
  end
  
  if verbose
%       fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik);
      fprintf(1, '.');
  end
  
  num_iter =  num_iter + 1;
  converged = em_converged(loglik, previous_loglik, thresh);
  previous_loglik = loglik;
%   LL = [LL loglik];

end

%%%%%%%%%
function [loglik, exp_num_trans, exp_num_visits1, lambda] = ...
    ess_phmm(prior, transmat, lambda_old, data)

% ESS_PHMM Compute the Expected Sufficient Statistics for a Poisson Hidden Markov Model.
%
% Outputs:
% exp_num_trans(i,j)   = sum_l sum_{t=2}^T Pr(Q(t-1) = i, Q(t) = j| Obs(l))
% exp_num_visits1(i)   = sum_l Pr(Q(1)=i | Obs(l))

% where Obs(l) = Obs(:,:,l) = O_1 .. O_T for sequence l
% Then 

verbose = 0;

numex = length(data);
O = size(data{1},1);
Q = length(prior);

exp_num_trans = zeros(Q,Q);
exp_num_visits1 = zeros(Q,1);
gamma_n = zeros(O,Q);
gammasum = zeros(Q,1);

loglik = 0;
if verbose, fprintf(1, 'forwards-backwards example # '); end
for ex=1:numex
  if verbose, fprintf(1, '%d ', ex); end
  %obs = data(:,:,ex);
  obs = data{ex};
  T = size(obs,2);
  
 
  B = poisson_prob(obs, lambda_old,1);
  
  
  [alpha, beta, gamma,  current_loglik, xi_summed] = fwdback(prior, transmat, B);
  
  
  loglik = loglik +  current_loglik; 
  if verbose, fprintf(1, 'll at ex %d = %f\n', ex, loglik); end

  exp_num_trans = exp_num_trans + xi_summed; % sum(xi,3);
  exp_num_visits1 = exp_num_visits1 + gamma(:,1);
  
  
  gamma_n = gamma_n + obs * gamma';
  gammasum = gammasum + sum(gamma, 2);

end
lambda = gamma_n ./ repmat(gammasum', [O 1]);

if verbose, fprintf(1, '\n'); end

