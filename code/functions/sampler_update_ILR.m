function [l_ast,b,rates,l_hood,rate_accept] = sampler_update_ILR(params,sample_old) 

rate_accept = repmat([0 realmin],params.n_params,1);

for r = 1:params.rep_MH  % Extra MH steps  
    n  = randperm(params.n_params);
    nn = n(1:params.n_joint);
    
    [sample_old.sigma_ast,sample_old.b,sample_old.rates,sample_old.l_hood,rate_accept] =...
                             MH_Joint_sampler(sample_old,params,nn,rate_accept);
                         
    combos = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];
    LL = zeros(1,params.n_params);
    LL(find(prod(sample_old.b(combos),2))) = 1;
    LL(unique(2*nchoosek(params.L,2)+find(sample_old.b))) = 1;
    LL(end) = 1;
    LL(nn)  = 0;
    
%     keyboard
    
    while sum(LL)>0              % All rates get updated
        % Choose order and pairings for updates
        ind=[];
        for l=1:params.n_joint
            % Randomly select a group of loads
            ind(l)=Discrete_sampler(LL);
            LL(ind(l))=0;
            if sum(LL)==0
                break;
            end
        end

        [sample_old.rates,sample_old.l_hood,rate_accept] =...
                     MH_sampler(sample_old,params,ind,rate_accept);
    end
    
    % Inactive rates
    LL = ones(1,params.n_params);
    LL(find(prod(sample_old.b(combos),2))) = 0;
    LL(unique(2*nchoosek(params.L,2)+find(sample_old.b))) = 0;
    LL(end) = 0;
    LL(nn)  = 0;
    INACT_rates = find(LL);

    prop_rates = gamrnd(params.phi,params.psi./params.phi,[params.n_params,1]); % from prior
    sample_old.rates(INACT_rates) = prop_rates(INACT_rates);
    rate_accept(INACT_rates,:) = rate_accept(INACT_rates,:) + 1;
    
end
l_ast  = sample_old.sigma_ast;
b      = sample_old.b;
rates  = sample_old.rates;
l_hood = sample_old.l_hood;

end

function [old_l_ast,old_b,old_rates,l_hood,rate_accept] = MH_Joint_sampler(sample_old,params,ind,rate_accept)

old_q      = sample_old.q;
old_l_ast  = sample_old.sigma_ast;
old_b      = sample_old.b;
old_rates  = sample_old.rates;

m    = sym(old_l_ast);
n    = sym(1:params.L);
old_kron = double(kroneckerDelta(n,m));
old_q_tilde = old_b.*params.xi + (1-old_b).*(1-params.xi);

prop_l_ast = Discrete_sampler(old_q./sum(old_q));

m    = sym(prop_l_ast);
n    = sym(1:params.L);
prop_kron = double(kroneckerDelta(n,m));

prop_b = binornd(ones(1,params.L),prop_kron + (1-prop_kron).*old_q_tilde);

prop_q_tilde =  prop_b.*params.xi + (1-prop_b).*(1-params.xi);

prop_rates = old_rates;
prop_rates(ind) = old_rates(ind)./params.p_lambda(ind).*randg(params.p_lambda(ind),[length(ind),1]);

prop_lhood             = sample_llhood(prop_l_ast,prop_b,prop_rates,params);
l_hood                 = sample_old.l_hood;

% keyboard 

b_priors = sum(prop_b.*log(prop_kron+(1-prop_kron).*old_q) ...
         + (1-prop_b).*log(1-(prop_kron+(1-prop_kron).*old_q)+realmin*prop_b)) ...
         - sum(old_b.*log(old_kron+(1-old_kron).*old_q) ...
         + (1-old_b).*log(1-(old_kron+(1-old_kron).*old_q)+realmin*old_b));
     
b_props  = sum(old_b.*log(old_kron+(1-old_kron).*prop_q_tilde) ...
         + (1-old_b).*log(1-(old_kron+(1-old_kron).*prop_q_tilde)+realmin*old_b)) ...
         - sum(prop_b.*log(prop_kron+(1-prop_kron).*old_q_tilde) ...
         + (1-prop_b).*log(1-(prop_kron+(1-prop_kron).*old_q_tilde)+realmin*prop_b));

rate_priors = sum((params.phi(ind)-1).*log((prop_rates(ind)+realmin)./(old_rates(ind)+realmin)) ...
            + (params.phi(ind)./params.psi(ind)).*(old_rates(ind) - prop_rates(ind)));      

rate_props  = sum((2*params.p_lambda(ind)-1).*(log(old_rates(ind)+realmin)-log(prop_rates(ind)+realmin)) ...
            + params.p_lambda(ind).*(prop_rates(ind)./(old_rates(ind)+realmin)-old_rates(ind)./(prop_rates(ind)+realmin)));

A = prop_lhood-l_hood+b_priors+b_props+rate_priors+rate_props;
    if  log(rand) < A*sample_old.omega_history
        old_l_ast             = prop_l_ast;
        old_b                 = prop_b;
        old_rates             = prop_rates;
        l_hood                = prop_lhood;
        rate_accept(ind,1)    = rate_accept(ind,1)+1 ;
    end
        rate_accept(ind,2)    = rate_accept(ind,2)+1 ;

end

function [old_rates,l_hood,rate_accept] = MH_sampler(sample_old,params,ind,rate_accept)

old_rates  = sample_old.rates;

prop_rates = old_rates;
prop_rates(ind) = old_rates(ind)./params.p_lambda(ind).*randg(params.p_lambda(ind),[length(ind),1]);

prop_lhood             = sample_llhood(sample_old.sigma_ast,sample_old.b,prop_rates,params);
l_hood                 = sample_old.l_hood;

% keyboard 

priors = (params.phi(ind)-1).*log((prop_rates(ind)+realmin)./(old_rates(ind)+realmin))+...
         (params.phi(ind)./params.psi(ind)).*(old_rates(ind) - prop_rates(ind));      

props  = (2*params.p_lambda(ind)-1).*(log(old_rates(ind)+realmin)-log(prop_rates(ind)+realmin))+...
         params.p_lambda(ind).*(prop_rates(ind)./(old_rates(ind)+realmin)-old_rates(ind)./(prop_rates(ind)+realmin));

A = prop_lhood-l_hood+sum(priors)+sum(props);
    if  log(rand) < A*sample_old.omega_history
        old_rates             = prop_rates;
        l_hood                = prop_lhood;
        rate_accept(ind,1)    = rate_accept(ind,1)+1 ;
    end
        rate_accept(ind,2)    = rate_accept(ind,2)+1 ;

end