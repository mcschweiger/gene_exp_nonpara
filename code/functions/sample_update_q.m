function [ q, rate_accept ] = sample_update_q(params,sample)
    rate_accept = repmat([0 realmin],params.L,1);

for r = 1:params.rep_q

    QQ       = ones(1,params.L);
    while sum(QQ) > 0
        q_upd = [];
        for qq = 1:params.q_joint
            q_upd(qq) = Discrete_sampler(QQ);
            QQ(q_upd(qq)) = 0;
            if sum(QQ) == 0
                break;
            end
        end
        [sample.q,rate_accept] =...
            MH_sampler(sample,params,q_upd,rate_accept);
    end
end
q      = sample.q;    

% keyboard
end

function [ q_old, rate_accept ] = MH_sampler(sample,params,q_upd,rate_accept)
    
    q_old     = sample.q;
    sigma_old = sample.sigma_ast;
    b_old     = sample.b;
    
    
    q_prop        = q_old;
    q_prop(q_upd) = betarnd(params.zeta/params.L,(params.L-1)/params.L,[1,length(q_upd)]);
    
    lp                  = q_upd;
    lp(lp == sigma_old) = [];
    
    m    = sym(sigma_old);
    n    = sym(q_upd);
    kron = double(kroneckerDelta(n,m));
% keyboard
    
    b_prior = b_old(q_upd).*(log(kron + (1-kron)*q_prop(q_upd))...
            - log(kron + (1-kron)*q_old(q_upd))) ...
            + (1-b_old(q_upd))*(log(1-(kron + (1-kron)*q_prop(q_upd))+eps*(b_old(q_upd))) ...
            - log(1-(kron + (1-kron)*q_old(q_upd))+eps*(b_old(q_upd))));

    sig_prior = log(q_prop(sigma_old)) - log(sum(q_prop)) ...
              - log(q_old(sigma_old))  + log(sum(q_old));

    logR  = sum(b_prior) + sig_prior;
    
%     keyboard

    if  log(rand) < logR*sample.omega_history
        q_old                     = q_prop;
        rate_accept(q_upd,1)      = rate_accept(q_upd,1)+1 ;
    end
    rate_accept(q_upd,2)    = rate_accept(q_upd,2)+1 ;

end