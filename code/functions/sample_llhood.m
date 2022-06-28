function [loglikelihood] = sample_llhood(sample_init,sample_b,sample_rates,params)

% keyboard

scale = params.t_f;

M = params.M;
g = sum(sample_b);

global A

% Chemical Master Equation Setup
A = double(generator_matrix(sample_b,sample_rates+log10(scale),params));

% ODE solver
tspan  = [params.t_0 params.t_f]./scale;

% keyboard
P_init = zeros(g*(M+1),1);

inits  = cumsum(sample_b);

P_init((M+1)*(inits(sample_init)-1)+1:(M+1)*(inits(sample_init)-1)+M+1) = params.Pm_init;

% OFF and zeros mRNA in nucleus and cytoplasm
options = odeset('RelTol',1e-4,'AbsTol',1e-5,'Jacobian',A,'JPattern',sparse(A),'MaxStep',1,'NonNegative',1);

% keyboard
clear persistent %<-- Reset persistent variables
try
    if all(diff(diff(params.obs_t)) < 1e-9)
        %tic
        times = params.obs_t./scale;
        P = zeros(length(P_init),length(params.obs_t));
        expo = expm(sparse((times(2)-times(1))*A));
        for i = 1:length(params.obs_t)
            t = times(i);
            if i == 1
                P(:,i) = P_init;
            else
                P(:,i) = expo*P(:,i-1);
            end
        end
        %t_exp = toc;
        
    else
        P  = ode15s(@dPdt, tspan, P_init, options);
        P  = deval(P,params.obs_t./scale);
    end
    %%
    % sum(P_g) per m value
    Psum = zeros(M+1,params.K);
    for i = 1:g
        Psum = Psum + P((i-1)*(M+1)+1:i*(M+1),:);
    end
    % keyboard
    %     if any(abs(1-sum(Psum)) > 1e-4)
    % %         keyboard
    %         disp('Warning: (1-P) > 1e-4')
    %     end
    %%
    % loglikelihood per cell and time
    loglikelihood_pertime = zeros(1,length(params.obs_t));
    
    for t = 1:length(params.obs_t)
        ind_loglikelihood = log(abs(Psum(params.r_ct{t}+1,t)));
        loglikelihood_pertime(t) = sum(ind_loglikelihood);
    end
    
    loglikelihood = sum(loglikelihood_pertime);
    % else
    %     loglikelihood = -inf;
    % end
catch
    %     keyboard
    disp(lasterr)
    loglikelihood = -inf;
end
clear A
end
