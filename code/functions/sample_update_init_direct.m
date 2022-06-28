function [ sigma_ast, b, l_hood] = sample_update_init_direct(params,sample)


for r = 1:params.rep_loads
% keyboard

    b = sample.b;
    q = sample.q;
    
logr = repmat(-inf,1,params.L);
for l = find(b)
    
    m    = sym(l);
    n    = sym(1:params.L);
    kron = double(kroneckerDelta(n,m));

    loglikelihood = sample_llhood(l,b,sample.rates,params);
    
    b_prior   = sum(b.*log(kron+(1-kron).*q)...
              + (1-b).*log(1-(kron+(1-kron).*q)+realmin*(b)));
    
    sig_prior = log(q(l)) - log(sum(q));

    logr(l)      = loglikelihood + b_prior + sig_prior;
              
end
% keyboard

sigma_ast    = Discrete_sampler( softmax(logr'*sample.omega_history));
b(sigma_ast) = 1;

l_hood = sample_llhood(sigma_ast,b,sample.rates,params);

end