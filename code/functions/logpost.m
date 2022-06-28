
function [logr] = logpost(l,params,L,BB,rates,q)
% Calculate the posterior of all possible loads
% keyboard

logr = nan(1,2^(L-1));
k    = 1;
if iscolumn(q)
    q = q';
end 
for bb = (l-1)*2^(L-1)+1:(l-1)*2^(L-1)+2^(L-1)

    loglikelihood = sample_llhood(l,BB(bb,:),rates,params);

    m    = sym(l);
    n    = sym(1:L);
    kron = double(kroneckerDelta(n,m));

    b_prior   = sum(BB(bb,:).*log(kron+(1-kron).*q)...
              + (1-BB(bb,:)).*log(1-(kron+(1-kron).*q)+eps*BB(bb,:)));

    l_prior = log(q(l)) - log(sum(q));

    logr(k)      = loglikelihood ...
                  + b_prior + l_prior;
              
    k = k+1;
end
end