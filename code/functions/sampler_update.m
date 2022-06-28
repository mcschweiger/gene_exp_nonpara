function [sample,rate_accept,logpost,tag] = sampler_update( sample, params )

%% update counter
sample.i = sample.i + 1;

%% Sample New ...


% keyboard
% Success probability vector
[sample.q, sample.accr_q]                          = sample_update_q(params,sample);

% Initial condition gene state and loads
[sample.sigma_ast, sample.b, sample.l_hood]        = sample_update_init_loads(params,sample);

% rates
[sample.rates,sample.l_hood,rate_accept,tag]       = sampler_update_rates(params,sample);
% % TESTING --------------
% rate_accept = repmat([0 realmin],params.n_params,1);
% tag = 1;
% % TESTING --------------

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

% keyboard

rate_priors = sum(-log(sqrt(2*pi)*params.psi(ACT_rates)).*...
         (-.5*(sample.rates(ACT_rates) - params.phi(ACT_rates)./params.psi(ACT_rates)).^2));
b_prior = sum(b(ll).*log(q(ll)) + (1-b(ll)).*log(1-q(ll)));
sigma_prior = log(q(sigma)) - log(sum(q));
q_prior     = gammaln((params.zeta-1)/params.L) - log(gamma(params.zeta/params.L)+gamma((params.L-1)/params.L))...
            + sum((params.zeta/params.L - 1).*log(q) + ((params.L-1)/params.L - 1).*(1-log(q)));

% keyboard

logpost = loglikelihood ...
        + b_prior + sigma_prior + q_prior + sum(rate_priors);


