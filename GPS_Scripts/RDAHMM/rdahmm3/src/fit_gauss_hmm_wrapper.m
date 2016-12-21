function [ll, pi, A, mu, sigma, Q] = fit_gauss_hmm_wrapper(data, N, params);
%
% Author: Robert Granat
%
% Wrapper function for fit_gauss_hmm, which is mex gateway to the RDAHMM
% code for training an HMM with Gaussian output distributions
%
% Inputs:
% 
% DATA: Time series of observations of length T and observation dimension D
% N: The number of hidden states in the HMM
% PARAMS [optional]: Parameters to use in the HMM training optimization 
% procedure.  Default parameters can be selected by submitting a string
% argument: 'vanilla' will perform standard Baum-Welsch EM with random
% initial states; 'kmeans' will perform Baum-Welsch but uses kmeans to
% select initial state parameters; 'regularize' will use Dirichlet and
% Euclidean distance priors to bias the EM solution (see Granat 2004); 
% 'anneal' will use regularized deterministic annealing EM to find the
% solution (see also Granat 2004)
%
% PARAMS can also be submitted as a structure with the following fields:
%
%           thresh [double]
%             peps [double]
%       regularize [0 or 1]
%         reg_type [0 or 1]
%        init_type [1 or 2]
%         omega_Q1 [double]
%         omega_Q2 [double]
%         omega_Q3 [double]
%      omega_sigma [double]
%           anneal [0 or 1]
%      anneal_step [double]
%    anneal_factor [double]
%         beta_min [double]
%         covgraph [1xN vector of ints]
%           ntries [int]
%         maxiters [int]
%             seed [int, select 0 for random initialization]
%
% For a full explanation of these parameters see the README for the 
% RDAHMM source code
%
% If PARAMS is not assigned, it defaults to 'vanilla'.
%
% Outputs:
%
% LL: log likelihood of the trained model
% PI: 1xN vector of initial model state probabilities
% A: NxN matrix of model state-to-state transition probabilities, A(i,j) is the
% probability of switching from state j to state i
% MU: DxN matrix of model output means, mu(:,i) is the outptut mean of the ith
% state
% SIGMA: Dx(D*N) matrix of model output covariances, sigma(:,(i-1)+1:3) is the
% covariance for the ith state
% Q: (1xT) vector of the individually most likely state assignments (integers)
%

if (nargin < 3)
  params = set_default_gauss_hmm_params(N);
end

if isstr(params)
  params = set_default_gauss_hmm_params(N, params);
end

[rows, cols] = size(data);

% Flip data matrix if it looks like it is oriented the wrong way
if rows >= cols 
  data = data';
end

[rows, cols] = size(data);

% Normalize data matrix values to lie between 0 and 1
min_data = min(data,[],2);
max_data = max(data,[],2);
range = max_data - min_data;
range_norm_data = (data - min_data * ones(1, cols)) ./ (range * ones(1,cols));

% Learn the HMM
[ll, pi, A, mu, sigma, Q] = fit_gauss_hmm(range_norm_data, N, params);

% Correct for the normalization of data values before returning
mu = mu .* (range * ones(1,N)) + min_data * ones(1,N);
for i = 0:(N-1)
  sigma(:,(i*rows+1):((i+1)*rows)) = ...
    sigma(:,(i*rows+1):((i+1)*rows)) .* (range * range'); 
end
