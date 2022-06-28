function   [params] = sampler_adapt_proposals(params,b,rates,i)

if isempty(rates)
    params.adapt_int        = 200;
    params.adapt_min        = params.adapt_int + 3e2;
    params.adapt_max        = inf;
elseif i > params.adapt_min && i < params.adapt_max
    [params] = sampler_adapt_shape(params,b,rates,i);
end

end

function  [params] = sampler_adapt_shape(params,b,rates,i)

params.cov = ( 2.5^2/params.n_params * cov(rates(:,i-params.adapt_int:i)') ...
    + 2.5^2/params.n_params * eps * eye(params.n_params));

if i == params.adapt_min+1
    disp('......................................................')
    disp('Adaptation of proposal covariance initiated')
    disp('......................................................')
    
elseif i == params.adapt_max-1
    disp('......................................................')
    disp('Adaptation of proposal covariance finished')
    disp('......................................................')
end
end



