function params = set_default_gauss_hmm_params(N, opt_type);

if (nargin < 2)
  opt_type = ['vanilla'];
end

if strcmp(opt_type, 'vanilla')
  params.thresh = 1.0e-3;
  params.peps = 1.0e-3;
  params.regularize = 0;
  params.reg_type = 0;
  params.init_type = 1;
  params.omega_Q1 = 0;
  params.omega_Q2 = 0;
  params.omega_Q3 = 0;
  params.omega_sigma = 1.0e-6;
  params.anneal = 0;
  params.anneal_step = 0;
  params.anneal_factor = 0;
  params.beta_min = 0;
  params.covgraph = zeros(1, N);
  params.ntries = 1;
  params.maxiters = 1000;
  params.seed = 0;
  params.viterbi = 0;
elseif strcmp(opt_type, 'kmeans')
  params.thresh = 1.0e-3;
  params.peps = 1.0e-3;
  params.regularize = 1;
  params.reg_type = 0;
  params.init_type = 2;
  params.omega_Q1 = 0;
  params.omega_Q2 = 0;
  params.omega_Q3 = 0;
  params.omega_sigma = 1.0e-6;
  params.anneal = 0;
  params.anneal_step = 0;
  params.anneal_factor = 0;
  params.beta_min = 0;
  params.covgraph = zeros(1, N);
  params.ntries = 10;
  params.maxiters = 1000;
  params.seed = 0;
  params.viterbi = 0;
elseif strcmp(opt_type, 'regularize')
  params.thresh = 1.0e-3;
  params.peps = 1.0e-3;
  params.regularize = 1;
  params.reg_type = 1;
  params.init_type = 1;
  params.omega_Q1 = 0.1;
  params.omega_Q2 = 0.1;
  params.omega_Q3 = 0.1;
  params.omega_sigma = 1.0e-6;
  params.anneal = 0;
  params.anneal_step = 0;
  params.anneal_factor = 0;
  params.beta_min = 0;
  params.covgraph = zeros(1, N);
  params.ntries = 10;
  params.maxiters = 1000;
  params.seed = 0;
  params.viterbi = 0;
elseif strcmp(opt_type, 'anneal')
  params.thresh = 1.0e-3;
  params.peps = 1.0e-3;
  params.regularize = 1;
  params.reg_type = 1;
  params.init_type = 1;
  params.omega_Q1 = 0;
  params.omega_Q2 = 0;
  params.omega_Q3 = 0.1;
  params.omega_sigma = 1.0e-6;
  params.anneal = 1;
  params.anneal_step = 0;
  params.anneal_factor = 1.1;
  params.beta_min = 0.1;
  params.covgraph = zeros(1, N);
  params.ntries = 10;
  params.maxiters = 1000;
  params.seed = 0;
  params.viterbi = 0;
end
