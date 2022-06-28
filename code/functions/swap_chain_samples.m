function [chain]= swap_chain_samples(chain,chain_idx,chain_order,tic_id,r)

% keyboard

a = max(chain_order(chain_idx), chain_order(chain_idx-1));
b = min(chain_order(chain_idx), chain_order(chain_idx-1));

delta_energy  = -chain(a).logpost(r)+ ...
                chain(b).logpost(r)                       ;
delta_omega   = chain(a).omega_history(r,:) - ...
                chain(b).omega_history(r,:)               ;

if ( delta_omega*delta_energy > log(rand(1)) )
% keyboard
disp(['Swap chains ' num2str(a) ' & ' num2str(b)])
% swap chains . . .
energy_tmp  = chain(b).logpost(r) ;
ll_tmp      = chain(b).loglike(r) ;
rates_tmp   = chain(b).rates(:,r);
b_tmp       = chain(b).loads(:,r);
l_tmp       = chain(b).sigma_ast(:,r);
q_tmp       = chain(b).q(:,r);

% Swap the necessary parameters

chain(b).logpost(r)     = chain(a).logpost(r) ;
chain(b).loglike(r)     = chain(a).loglike(r) ;
chain(b).rates(:,r)     = chain(a).rates(:,r);
chain(b).loads(:,r)     = chain(a).loads(:,r);
chain(b).sigma_ast(:,r) = chain(a).sigma_ast(:,r);
chain(b).q(:,r)         = chain(a).q(:,r);

chain(a).logpost(r)     = energy_tmp;
chain(a).loglike(r)     = ll_tmp;
chain(a).rates(:,r)     = rates_tmp;
chain(a).loads(:,r)     = b_tmp;
chain(a).sigma_ast(:,r) = l_tmp;
chain(a).q(:,r)         = q_tmp;

% Update the samples
smp_energy_tmp  = chain(b).sample.logpost ;
smp_ll_tmp      = chain(b).sample.l_hood ;
smp_rates_tmp   = chain(b).sample.rates;
smp_b_tmp       = chain(b).sample.b;
smp_l_tmp       = chain(b).sample.sigma_ast;
smp_q_tmp       = chain(b).sample.q;

chain(b).sample.logpost   = chain(a).sample.logpost ;
chain(b).sample.l_hood    = chain(a).sample.l_hood ;
chain(b).sample.rates     = chain(a).sample.rates;
chain(b).sample.b         = chain(a).sample.b;
chain(b).sample.sigma_ast = chain(a).sample.sigma_ast;
chain(b).sample.q         = chain(a).sample.q;

chain(a).sample.logpost   = smp_energy_tmp ;
chain(a).sample.l_hood    = smp_ll_tmp ;
chain(a).sample.rates     = smp_rates_tmp;
chain(a).sample.b         = smp_b_tmp;
chain(a).sample.sigma_ast = smp_l_tmp;
chain(a).sample.q         = smp_q_tmp;

% keyboard



% Record acceptance
chain(b).swap_acceptance(r,:)   = 1;
end

chain(a).swap_time(1,:) = toc(tic_id);  