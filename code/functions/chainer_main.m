function [chain]= chainer_main(chain_init,d_length,nchains,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chainer_main will manipulate the MCMC chain at hand.
%
%
% to initialize (create a new chain), all you need is options:
% chain = chainer_main([]   ,  0, opts, true, []  );
% to expand (generate MCMC samples), all you need is a chain:
% chain = chainer_main(chain,+25, []  , true, true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')
flag_status = 1;

% Initialize a chain
if d_length == 0
    tic_id = tic;
    % MCMC
    chain(nchains) = struct();
    parfor chain_idx = 1:nchains
        
        chain(chain_idx).params = chainer_init_params(opts);%set up params
        chain(chain_idx).params = sampler_adapt_proposals(chain(chain_idx).params,[],[],[]);
        
        chain(chain_idx).length = 1;                        %monitor chain length
        
        chain(chain_idx).sizeGB = nan;                      %monitor chain memory size
        
        chain(chain_idx).sample = [];                       %store current sample
        chain(chain_idx).sample = chainer_init_sample(chain_idx,chain(chain_idx).params,opts);
        
        
        
        
        % Store the history of the chain
        chain(chain_idx).i         = cast( chain(chain_idx).sample.i    , 'uint64') ;
        
        chain(chain_idx).q         = cast( chain(chain_idx).sample.q', 'single') ;
        
        chain(chain_idx).loads     = cast( chain(chain_idx).sample.b', 'single') ;
        
        chain(chain_idx).sigma_ast = cast( chain(chain_idx).sample.sigma_ast, 'single') ;
        
        chain(chain_idx).rates     = cast( chain(chain_idx).sample.rates, 'single') ;
        
        chain(chain_idx).loglike   = cast( chain(chain_idx).sample.l_hood, 'single') ;
        
        chain(chain_idx).logpost   = cast( chain(chain_idx).sample.logpost, 'single') ;
        
        chain(chain_idx).accr_rates =  cast( chain(chain_idx).sample.accr_rates    , 'single') ;
        
        chain(chain_idx).swap_acceptance =  cast( chain(chain_idx).sample.swap_acceptance    , 'single') ;
        
        chain(chain_idx).relstep_history =  cast( chain(chain_idx).sample.relstep_history'    , 'single') ;
        
        chain(chain_idx).omega_history   =  cast( chain(chain_idx).sample.omega_history      , 'single') ;
        
        chain(chain_idx).swap_time       =  cast( chain(chain_idx).sample.swap_time          , 'single') ;
        
        %display some initialization info
        if flag_status
            disp(['CHAINER for ', 'chain_id= ',num2str(chain_idx), ' initiated (total time = ',num2str(toc(tic_id)),' s)'])
        end
        
    end
    
    chain_order   = randperm(nchains);      %propose random chain switches
    for chain_idx = 2:length(chain_order)   %then swap them, if accepted
        [chain]= swap_chain_samples(chain,chain_idx,chain_order,tic_id,1);
    end
    
% Expand a chain
elseif d_length > 0
    tic_id = tic;
    chain(nchains) = struct();
    parfor chain_idx = 1:nchains

        
        chain(chain_idx).params = chain_init(chain_idx).params;
        chain(chain_idx).length = chain_init(chain_idx).length + d_length;

        
        chain(chain_idx).sample    = chain_init(chain_idx).sample;
        chain(chain_idx).i         = [ chain_init(chain_idx).i      zeros( 1,d_length       , 'like', chain_init(chain_idx).i)];
        chain(chain_idx).q         = [ chain_init(chain_idx).q     nan( chain(chain_idx).params.L, d_length,'like',chain_init(chain_idx).q  )];
        chain(chain_idx).loads     = [ chain_init(chain_idx).loads nan( chain(chain_idx).params.L, d_length,'like',chain_init(chain_idx).loads  )];
        chain(chain_idx).sigma_ast = [ chain_init(chain_idx).sigma_ast  zeros( 1,d_length       , 'like', chain_init(chain_idx).sigma_ast)];
        chain(chain_idx).rates     = [ chain_init(chain_idx).rates nan( chain(chain_idx).params.n_params, d_length,'like',chain_init(chain_idx).rates  )];
        chain(chain_idx).logpost   = [ chain_init(chain_idx).logpost nan( 1, d_length,'like',chain_init(chain_idx).logpost  )];
        chain(chain_idx).loglike   = [ chain_init(chain_idx).loglike nan( 1, d_length,'like',chain_init(chain_idx).loglike  )];
        
        chain(chain_idx).accr_rates =  [ chain_init(chain_idx).accr_rates ];
        
        chain(chain_idx).swap_acceptance =  [ chain_init(chain_idx).swap_acceptance    ; zeros( d_length,1 ,'like', chain_init(chain_idx).swap_acceptance)];
        chain(chain_idx).relstep_history =  [ chain_init(chain_idx).relstep_history   ; zeros( d_length,chain_init(chain_idx).params.nparams ,'like', chain_init(chain_idx).relstep_history)];
        chain(chain_idx).omega_history   =  [ chain_init(chain_idx).omega_history      ; zeros( d_length,1 ,'like', chain_init(chain_idx).omega_history)];
        chain(chain_idx).swap_time       =  [ chain_init(chain_idx).swap_time          ; zeros( d_length,1 ,'like', chain_init(chain_idx).swap_time)];
        

    end

    
    for chain_idx = 1:nchains
        r = chain_init(chain_idx).length+1;
    end
    
    while r <= chain(chain_idx).length

        parfor chain_idx = 1 : nchains

            [chain(chain_idx).sample,chain(chain_idx).sample.accr_rates,chain(chain_idx).sample.logpost, chain(chain_idx).tag ] = ...
                sampler_update(chain(chain_idx).sample,chain(chain_idx).params);
            
            chain(chain_idx).accr_rates = chain(chain_idx).accr_rates + chain(chain_idx).sample.accr_rates;
            
            
            
            
            if mod(chain(chain_idx).sample.i,chain(chain_idx).params.i_skip) == 0
                
                chain(chain_idx).i(r)            = chain(chain_idx).sample.i      ;
                chain(chain_idx).q(:,r)          = chain(chain_idx).sample.q;
                chain(chain_idx).loads(:,r)      = chain(chain_idx).sample.b;
                chain(chain_idx).sigma_ast(r)    = chain(chain_idx).sample.sigma_ast;
                chain(chain_idx).rates(:,r)      = chain(chain_idx).sample.rates  ;
                chain(chain_idx).logpost(r)      = chain(chain_idx).sample.logpost;
                chain(chain_idx).loglike(r)      = chain(chain_idx).sample.l_hood;
                
                chain(chain_idx).params          = ...
                    sampler_adapt_proposals(chain(chain_idx).params,...
                    chain(chain_idx).loads,...
                    chain(chain_idx).rates,...
                    chain(chain_idx).sample.i);
                
                
                chain(chain_idx).accr_rates = chain(chain_idx).sample.accr_rates;
                
                chain(chain_idx).swap_acceptance(r,:) = chain(chain_idx).sample.swap_acceptance  ;
                chain(chain_idx).relstep_history(r,:) = chain(chain_idx).params.p_lambda;
                chain(chain_idx).omega_history(r,:)   = chain(chain_idx).sample.omega_history     ;
                chain(chain_idx).swap_time(r,:)       = chain(chain_idx).sample.swap_time        ;
                
                
                
            end
        end
        %% Try to swap chains        
        if mod(chain(chain_idx).sample.i,chain(chain_idx).params.chain_MCMCsteps) == 0
            chain_order   = randperm(nchains);
            for chain_idx = 2:2:length(chain_order)
                [chain]= swap_chain_samples(chain,chain_idx,chain_order,tic_id,r);
            end
        end
        
        

        for chain_idx = 1 : nchains
            %every 10 samples, display some info about the chain
            if flag_status && mod(chain(chain_idx).sample.i,10) == 0
                

                combos = [nchoosek(1:chain(chain_idx).params.L,2); flip(nchoosek(1:chain(chain_idx).params.L,2),2)];
                LL                                                = zeros(1,chain(chain_idx).params.n_params);
                LL(find(prod(chain(chain_idx).sample.b(combos),2)))                = 1;
                LL(unique(2*nchoosek(chain(chain_idx).params.L,2)+find(chain(chain_idx).sample.b))) = 1;
                LL(end)                                           = 1;
                ACT_rates                                         = find(LL);
                
                s_type = [{'MH S'},{'MH J'},{'HMC'},{'IBR J'}]; 
                
                
                disp([  'i = ', num2str(chain(chain_idx).sample.i,'%d'),...
                    ' - ',s_type{chain(chain_idx).tag},' acc_rate= ' ,'(',...
                    num2str(chain(chain_idx).accr_rates(ACT_rates,1)'./...
                    chain(chain_idx).accr_rates(ACT_rates,2)' * 100, ' %#6.2f'),...
                    ')' ]);
                disp(['b = ' num2str(chain(chain_idx).sample.b)])
            end
            
            
            if (mod(r, chain(chain_idx).params.chain_MCMCsteps)==0 &&  r <= chain(chain_idx).params.adapt_last )
                [chain(chain_idx).sample.omega_history] = sampler_adapt_omega( chain(chain_idx).omega_history(r,:),chain, r, chain(chain_idx).params, chain_idx );
                
            end
        end
        
        r = r+1;
    end
end



end

%% auxiliary functions

function sizeGB = get_sizeGB(chain)
sizeGB = whos( inputname(1) );
sizeGB = sizeGB.bytes/1024^3;
end


% ---- adapt max temperature -----------------------------------------

function [omega] = sampler_adapt_omega( omega,chain, swap_idx, opts, chain_idx )
% adapt omega
fprintf(1,'------------------------------------------------------\n');
fprintf(1,' Adapting chain temperatures . . .\n');
% keyboard
% based on acceptance rate of all chains
old_omega = omega;
% chain 1 has fixed temperature. adjust chains in order of coolest to hottest.
%   whenever a temperature changes, we proportionally increase the temperature of the hotter chains as well
if chain_idx > 1
    % compute adaption factor
    %     adaption_factor = ((sum( chain(chain_idx-1).swap_acceptance((swap_idx-opts.adapt_omega_interval+1):swap_idx) ) ...
    %                         / floor(swap_idx/opts.chain_MCMCsteps))/opts.optimal_swap_acceptance)^opts.adapt_omega_rate;
    adaption_factor = ((sum( chain(chain_idx-1).swap_acceptance(1:swap_idx) ) ...
        / floor(swap_idx/opts.chain_MCMCsteps))/opts.optimal_swap_acceptance)^opts.adapt_omega_rate;
    adaption_factor = min( [max([opts.min_adaptation_factor, adaption_factor]), opts.max_adaptation_factor] );
    % make sure omega doesn't get too close to its cooler neighbor
    % keyboard
    adaption_factor = max( [adaption_factor, 2*omega/(omega + chain(chain_idx-1).omega_history(swap_idx,:))] );
    % multiply this chain and higher temperature chains by adaption_factor
    
    omega = omega/adaption_factor;
    
    fprintf(1,'  chain %d: recent acceptance=%-8.3g old temperature=%-6.4f new temperature=%-6.4f\n', ...
        chain_idx, (sum( chain(chain_idx-1).swap_acceptance((swap_idx-opts.adapt_omega_interval+1):swap_idx) ) ...
        / opts.adapt_omega_interval), 1/old_omega, 1/omega );
end
end





