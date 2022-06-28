function [sigma_ast, b, l_hood] = sample_update_init_loads(params,sample)
% keyboard
b     = sample.b;
q     = sample.q;
rates = sample.rates;

L = params.L;

% keyboard

b_tmp = ones(1,L);
BB    =  repmat(b_tmp,L*2^(L-1),1);
for l = 1:params.L
    lp = 1:L;
    lp(l) = [];
    
    lp_size =length(lp);

    BB((l-1)*2^(L-1)+1:(l-1)*2^(L-1)+2^(L-1),lp) = load_sets(lp_size); 
end


logp(1:L)  = struct();
for l = 1:params.L        
    logp(l).lpost = logpost(l,params,L,BB,rates,q);
end

logr = [];
for l = 1:L
    logr = [logr logp(l).lpost];
end
    
%     keyboard
for r = 1:params.rep_loads   
    load_idx = Discrete_sampler( softmax(logr'*sample.omega_history));
    l_idx    = find(load_idx >= (1:L-1)*2^(L-1)+1);
    if isempty(l_idx)
        l_idx = 0;
    end
    
    sigma_ast = l_idx(end)+1;
    

    
    b = BB(load_idx,:);
    l_hood = sample_llhood(sigma_ast,b,rates,params);
    
end
end
