function [ phi_ast, rate_accept ] = sample_update_phi(params,sample)
    rate_accept = [0 realmin];

    for r = 1:params.rep_phi
        [ sample.phi_ast, rate_accept ] = MH_sampler(sample,params,rate_accept);
    end
    
    phi_ast = sample.phi_ast ;
end

function [ phi_old, rate_accept ] = MH_sampler(sample,params,rate_accept)
%     keyboard

    phi_old   = sample.phi_ast;
    b_old     = sample.b;
    sigma_old = sample.sigma_ast;
    
    phi_prop  = dirrnd(params.alpha*phi_old + params.kappa);
    
    logR      = log(b_old(sigma_old)*phi_prop(sigma_old)) - log(sum(b_old.*phi_prop)) ...
              - log(b_old(sigma_old)*phi_old(sigma_old)) - log(sum(b_old.*phi_old)) ...
              + sum((params.alpha*params.psi_ast' - 1).*log(phi_prop+realmin)) ...
              + sum((params.alpha*params.psi_ast' - 1).*log(phi_old+realmin)) ...
              - log(prod(gamma(params.alpha*phi_prop+params.kappa))) + sum((params.alpha*phi_prop - 1).*log(phi_old+realmin)) ...
              + log(prod(gamma(params.alpha*phi_old+params.kappa))) - sum((params.alpha*phi_old - 1).*log(phi_prop+realmin));
    if  log(rand) < logR
        phi_old             = phi_prop;
        rate_accept(1)    = rate_accept(1)+1 ;
    end
        rate_accept(2)    = rate_accept(2)+1 ;
% keyboard
end