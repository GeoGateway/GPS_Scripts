function [ll, Q] = label_using_gauss_hmm_wrapper(data, pi, A, mu, sigma, range, min_data, viterbi, addstate);

[rows, cols] = size(data);

% Flip data matrix if it looks like it is oriented the wrong way
if rows >= cols
  data = data';
end

[rows, cols] = size(data);

if nargin <= 5
  min_data = min(data,[],2);
  max_data = max(data,[],2);
  range = max_data - min_data;
end

% Normalize data matrix values in accordance with training data 
range_norm_data = (data - min_data * ones(1, cols)) ./ (range * ones(1,cols));

N = length(pi);

% Correct for the normalization of data 
mu = (mu - min_data * ones(1,N)) ./ (range * ones(1,N));
for i = 0:(N-1)
  sigma(:,(i*rows+1):((i+1)*rows)) = ...
    sigma(:,(i*rows+1):((i+1)*rows)) ./ (range * range');
end

if (nargin <= 7)
  hmm.N = length(pi);
  hmm.D = size(mu,1);
  hmm.pi = pi;
  hmm.A = A;
  hmm.mu = mu;
  hmm.sigma = sigma;
  viterbi = 0;
elseif (nargin <= 9)
  if addstate
    hmm.N = N + 1;
    hmm.D = size(mu,1);
    hmm.pi = [(N * pi / (N+1)) 1/(N+1)];
    hmm.A = ones(N+1);
    hmm.A(1:N,1:N) = N * A;
    hmm.A = hmm.A / (N+1);
    hmm.mu = [mu (mean(data,2) ./ range + min_data)];
    sigma = [sigma cov(data') ./ (range * range')];
  end
end  

[ll, Q] = label_using_gauss_hmm(range_norm_data, hmm, viterbi);
