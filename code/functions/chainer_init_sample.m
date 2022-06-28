function sample = chainer_init_sample(chain_idx,params,opts)
% Compute initial sample and its likelihood

sample.i  = 0;

% Initial proposals--------------------------------------------------------
sample.params = opts;

%Sample the model and its parameters from the specified priors
sample.q         = betarnd(params.zeta/params.L,(params.L-1)/params.L,[1,params.L]);
sample.sigma_ast = Discrete_sampler(sample.q./sum(sample.q));
sample.b(sample.sigma_ast) = 1;
lp                         = 1:params.L;
lp(lp == sample.sigma_ast) = [];
sample.b(lp)               = binornd(ones(1,params.L-1),sample.q(lp));
sample.rates = normrnd(params.phi,params.psi,[params.n_params,1]);


%% Readouts for monitoring
disp('...................................................')
disp(['Initial gene state: ', num2str(sample.sigma_ast)])
disp('...................................................')
disp(['Initial loads: ', num2str(sample.b)])
disp('...................................................')
for i = 1:params.n_params
    disp(['Initial ',char(params.names(i)),...
          blanks(max(length(char(params.names))) - length(char(params.names(i)))),...
          ' : ', num2str(sample.rates(i))])
end

disp('...................................................')

%Trivial acceptance rates
sample.accr_q     = repmat([0 realmin],params.L,1);
sample.accr_sigma = repmat([0 realmin],params.L,1);
sample.accr_loads = repmat([0 realmin],params.L,1);
sample.accr_rates = repmat([0 realmin],params.n_params,1);


%Sample log likelihood
[sample.l_hood] = sample_llhood(sample.sigma_ast,sample.b,sample.rates,params);

sample.swap_acceptance = 0;                                                    % number of accepted swaps (per chain per swap)
sample.omega_history   = params.max_omega * (params.omega_init)^(chain_idx-1); % omega = inverse chain temperature
sample.relstep_history = params.relstep_init./sample.omega_history ;           % relative step size w.r.t. log-parameter interval
sample.swap_time       = 0;                                                    % time per swap

[sample.l_hood] = sample_llhood(sample.sigma_ast,sample.b,sample.rates,params);

%% Calculate the log posterior

loglikelihood = sample.l_hood;

combos = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];
LL                                                = zeros(1,params.n_params);
LL(find(prod(sample.b(combos),2)))                = 1;
LL(unique(2*nchoosek(params.L,2)+find(sample.b))) = 1;
LL(end)                                           = 1;
ACT_rates                                         = find(LL);

ll = 1:params.L;
ll(ll == sample.sigma_ast) = [];

q     = sample.q;
b     = sample.b;
sigma = sample.sigma_ast;


rate_priors = sum(-log(sqrt(2*pi)*params.psi(ACT_rates)).*...
         (-.5*(sample.rates(ACT_rates) - params.phi(ACT_rates)./params.psi(ACT_rates)).^2));
b_prior = sum(b(ll).*log(q(ll)) + (1-b(ll)).*log(1-q(ll)));
sigma_prior = log(q(sigma)) - log(sum(q));
q_prior     = gammaln((params.zeta-1)/params.L) - log(gamma(params.zeta/params.L)+gamma((params.L-1)/params.L))...
            + sum((params.zeta/params.L - 1).*log(q) + ((params.L-1)/params.L - 1).*(1-log(q)));


sample.logpost = loglikelihood ...
        + b_prior + sigma_prior + q_prior + sum(rate_priors);


end

