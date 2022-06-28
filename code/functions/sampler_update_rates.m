function [rates,l_hood,rate_accept,tag] = sampler_update_rates(params,sample_old)

rate_accept = repmat([0 realmin],params.n_params,1);

                                       % Randomly choose single or joint sampling
rep = params.rep;


if sample_old.i < params.adapt_min
    tag = randi(2);
else
    tag = randi(3);
end


if sample_old.i < params.adapt_max
    rep = 1;
end

% keyboard
switch tag
    case 1
% Single-------------------------------------------------------------------
        combos = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];
        LL = zeros(1,params.n_params);
        LL(find(prod(sample_old.b(combos),2))) = 1;
        LL(unique(2*nchoosek(params.L,2)+find(sample_old.b))) = 1;
        LL(end) = 1;
        ACT_rates = find(LL);
%         ACT_rates = [24];
        INACT_rates = find(LL == 0);
        for r = 1:rep  % Extra samples
            % MH sampling for active rates
            update_order = randperm(size(ACT_rates,2));
            for i = 1:size(ACT_rates,2)
                ind  = ACT_rates(update_order(i));
                [sample_old.rates,sample_old.l_hood,rate_accept] =...
                        MH_sampler(sample_old,params,ind,rate_accept);
            end
            % Direct sampling for inactive rates
            prop_rates = normrnd(params.phi,params.psi,[params.n_params,1]); % from prior
            sample_old.rates(INACT_rates) = prop_rates(INACT_rates);
            rate_accept(INACT_rates,:) = rate_accept(INACT_rates,:) + 1;

        end

    case 2                           % Joint sampling after adapting
% Joint--------------------------------------------------------------------
       for r = 1:rep  % Extra MH steps
            % Active rates
            combos = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];
            LL = zeros(1,params.n_params);
            LL(find(prod(sample_old.b(combos),2))) = 1;
            LL(unique(2*nchoosek(params.L,2)+find(sample_old.b))) = 1;
            LL(end) = 1;

            ind = find(LL);
%             ind = [24];

            [sample_old.rates,sample_old.l_hood,rate_accept] =...
                         MH_sampler(sample_old,params,ind,rate_accept);
            % Inactive rates
            LL = ones(1,params.n_params);
            LL(find(prod(sample_old.b(combos),2))) = 0;
            LL(unique(2*nchoosek(params.L,2)+find(sample_old.b))) = 0;
            LL(end) = 0;
            INACT_rates = find(LL);

            prop_rates = normrnd(params.phi,params.psi,[params.n_params,1]); % from prior
            sample_old.rates(INACT_rates) = prop_rates(INACT_rates);
            rate_accept(INACT_rates,:) = rate_accept(INACT_rates,:) + 1;



       end
    case 3
        for r = 1:rep  % Extra MH steps
            % Active rates
            combos = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];
            LL = zeros(1,params.n_params);
            LL(find(prod(sample_old.b(combos),2))) = 1;
            LL(unique(2*nchoosek(params.L,2)+find(sample_old.b))) = 1;
            LL(end) = 1;
            ACT_rates = find(LL);

%             keyboard

            [sample_old.rates,sample_old.l_hood,rate_accept] =...
                         sampler_HMC(sample_old,params,ACT_rates,rate_accept);

            % Inactive rates
            LL = ones(1,params.n_params);
            LL(find(prod(sample_old.b(combos),2))) = 0;
            LL(unique(2*nchoosek(params.L,2)+find(sample_old.b))) = 0;
            LL(end) = 0;
            INACT_rates = find(LL);

            prop_rates = normrnd(params.phi,params.psi,[params.n_params,1]); % from prior
            sample_old.rates(INACT_rates) = prop_rates(INACT_rates);
            rate_accept(INACT_rates,:) = rate_accept(INACT_rates,:) + 1;

       end
end
rates  = sample_old.rates;
l_hood = sample_old.l_hood;
end

function [old_rates,l_hood,rate_accept] = MH_sampler(sample_old,params,ind,rate_accept)

old_rates  = sample_old.rates;

prop_rates = old_rates;
prop_rates(ind) = mvnrnd(old_rates(ind)',params.cov(ind,ind));

prop_lhood             = sample_llhood(sample_old.sigma_ast,sample_old.b,prop_rates,params);
l_hood                 = sample_old.l_hood;

lpriors   = .5*(((old_rates(ind) - params.phi(ind)).^2 - (prop_rates(ind) - params.phi(ind)).^2)./params.psi(ind).^2);

% lprops    = -.5*(old_rates(ind)-prop_rates(ind))'*inv(params.cov(ind,ind))*(old_rates(ind)-prop_rates(ind))...
%           + .5*(prop_rates(ind)-old_rates(ind))'*inv(params.cov(ind,ind))*(prop_rates(ind)-old_rates(ind));

A = prop_lhood-l_hood+sum(lpriors);

% keyboard
    if  log(rand) < A*sample_old.omega_history
        old_rates             = prop_rates;
        l_hood                = prop_lhood;
        rate_accept(ind,1)    = rate_accept(ind,1)+1 ;
    end
        rate_accept(ind,2)    = rate_accept(ind,2)+1 ;

end
