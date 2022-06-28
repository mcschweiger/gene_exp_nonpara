function [rates,lhood,rate_accept] = sampler_HMC(sample,params,ACT_rates,rate_accept)

rates        = sample.rates;
omega        = sample.omega_history;

rates_old = rates(ACT_rates);

nparams   = size(ACT_rates,2);

h        = params.HMC_eps*randg();
HMC_L    = params.HMC_L;
demo     = false;

% keyboard
m_MSp        = gamrnd(params.phi_m',params.psi_m'./params.phi_m',[1,nparams]);

% MOMENTUM
sample_MSp  = m_MSp.*randn(1,nparams);

% Strang-Splitting PREP
propos_mSq  = nan(HMC_L,nparams);
propos_MSp  = nan(HMC_L,nparams);

% FIRST STEP
i = 1;
propos_mSq(i,:)   = rates_old;
propos_MSp(i,:)   = sample_MSp;

[V_q]            = find_V_grad(rates_old,sample,params);
V_q              = V_q .* (10.^rates_old'*log(10)) * omega;

% keyboard
    if all(~isnan(V_q))

        if demo;Gim = HMC_visual([],nparams,HMC_L,propos_mSq,propos_MSp); end

        for i = 2:HMC_L

            % Half-step (Applying only V_q)
            propos_MSp(i,:) = propos_MSp(i-1,:) - 0.5 * h * omega * V_q;

            % Whole-step (Applying M and L_q)
            propos_mSq(i,:) = ((4*h*omega.*params.psi(ACT_rates)'.^2.*(1./m_MSp)).*propos_MSp(i,:)...
                                -(2*h^2*omega^2*params.phi(ACT_rates)').*(1./m_MSp) + (4*params.psi(ACT_rates)'.^2 + h^2*omega^2*(1./m_MSp)).*propos_mSq(i-1,:))./...
                                (4*params.psi(ACT_rates)'.^2 -  h^2*omega^2*(1./m_MSp));

            propos_MSp(i,:) = propos_MSp(i,:) + h*omega*(((propos_mSq(i-1,:)+propos_mSq(i,:))./2 - params.phi(ACT_rates)')./(params.psi(ACT_rates).^2)');

            [V_q]            = find_V_grad(propos_mSq(i,:),sample,params);
            V_q = V_q .* (10.^propos_mSq(i,:)*log(10));
            
            if isnan(V_q)
               break
            end

            % Half-step (Applying only V_q)
            propos_MSp(i,:) = propos_MSp(i,:) - 0.5 * h * omega * V_q;

            if demo;HMC_visual(Gim,nparams,HMC_L,propos_mSq,propos_MSp); end
        end
    end

    if all(~isnan(V_q))
        L_old      = sum(.5 * ((propos_mSq(1,:)- params.phi(ACT_rates)')./params.psi(ACT_rates)').^2);
        T_old      =  sum(0.5* (sample_MSp.^2)./m_MSp);
        H_old      = omega*(find_H1(propos_mSq(1,:),ACT_rates,sample,params) + L_old + T_old);

        L_prop      = sum(.5 * ((propos_mSq(end,:)- params.phi(ACT_rates)')./params.psi(ACT_rates)').^2);
        T_prop      = sum(0.5* (propos_MSp(end,:).^2)./m_MSp);
        H_prop      = omega*(find_H1(propos_mSq(end,:),ACT_rates,sample,params) + L_prop + T_prop);


        log_a = H_old - H_prop;
    else
        log_a = -inf;
    end

        if log(rand) < log_a
            rates            = sample.rates;
            rates(ACT_rates) = propos_mSq(end,:)';
            lhood            = sample_llhood(sample.sigma_ast,sample.b,rates,params);
            rate_accept(ACT_rates,1) = rate_accept(ACT_rates,1)+1;
        else
            rates = sample.rates;
            lhood = sample_llhood(sample.sigma_ast,sample.b,rates,params);
        end
            rate_accept(ACT_rates,2) = rate_accept(ACT_rates,2)+1;

        if ~iscolumn(rates_old)
        rates_old = rates_old';
        end
end

%% plotter
function Gim = HMC_visual(Gim,n_params,HMC_L,propos_mS_q,...
                                     propos_MS_p)

if isempty(Gim)
    figure(88)

    col = [0 1 0;%green
           0 0 1;%blue
           1 0 0];%red

    subplot(1,2,1)
    Gim.ax_Di{1} = plot(1:HMC_L,propos_mS_q,'o-'        ,'color',col(1,:));
    xlim([0 HMC_L+1])

    subplot(1,2,2)
    Gim.ax_Dm{1} = plot(propos_MS_p,propos_mS_q,'o-'       ,'color',col(1,:));

end

for n=1:n_params
    Gim.ax_Di{1}(n).YData = propos_mS_q(:,n);
    Gim.ax_Dm{1}(n).XData = propos_MS_p(:,n);
    Gim.ax_Dm{1}(n).YData = propos_mS_q(:,n);
end

drawnow

end


%% HMC potential
function V = find_H1(rates,ACT_rates,sample,params)

M = params.M;
g = sum(sample.b);

scale = params.t_f;

sample_rates = repmat(-1,params.n_params,1);
sample_rates(ACT_rates) = rates;

global A

% Chemical Master Equation Setup
A = generator_matrix(sample.b,sample_rates+log10(scale),params);

% ODE solver
tspan  = [params.t_0 params.t_f]./scale;
% keyboard
P_init = zeros(g*(M+1),1);

inits  = cumsum(sample.b);

P_init((M+1)*(inits(sample.sigma_ast)-1)+1:(M+1)*(inits(sample.sigma_ast)-1)+M+1) = params.Pm_init;

% OFF and zeros mRNA in nucleus and cytoplasm
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'Jacobian',A,'JPattern',sparse(A),'MaxStep',1,'NonNegative',1);

% keyboard
clear persistent %<-- Reset persistent variables
try
    P  = ode15s(@dPdt, tspan, P_init, options);
    P  = deval(P,params.obs_t./scale);

    %%
    % sum(P_g) per m value
    Psum = zeros(M+1,params.K);
    for i = 1:g
        Psum = Psum + P((i-1)*(M+1)+1:i*(M+1),:);
    end
% keyboard
%     if any(abs(1-sum(Psum)) > 1e-4)
%         disp('Warning: (1-P) > 1e-4')
%     end
    %%
    % loglikelihood per cell and time
    loglikelihood_pertime = zeros(1,length(params.obs_t));

    for t = 1:length(params.obs_t)
        ind_loglikelihood = log(abs(Psum(params.r_ct{t}+1,t)));
        loglikelihood_pertime(t) = sum(ind_loglikelihood);
    end

    V = -sum(loglikelihood_pertime);
catch
    disp('ode15s warning ignored, sample rejected')
    V = nan;
end
clear A
end


%% HMC potential gradient
function [V_mS] = find_V_grad(rates,sample,params)
if iscolumn(rates)
    rates = rates';
end

global M b L

scale = params.t_f;
% scale = 1;

M = params.M;
b = sample.b;

g = sum(sample.b);

P_init = zeros(g*(M+1),1);

inits  = cumsum(sample.b);

P_init((M+1)*(inits(sample.sigma_ast)-1)+1:(M+1)*(inits(sample.sigma_ast)-1)+M+1) = params.Pm_init;

% keyboard
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
try
    [~,P,dPdr] = sens_sys('CME',params.obs_t./scale,P_init,opts,rates'+log10(scale));


    % keyboard

    Psum    = zeros(params.K,M+1);
    dPdrsum = zeros([size(dPdr,1),M+1,size(dPdr,3)]);
    for i = 1:g
        Psum      = Psum + P(:,(i-1)*(M+1)+1:i*(M+1));
        dPdrsum   = dPdrsum + dPdr(:,(i-1)*(M+1)+1:i*(M+1),:);
    end

    sumj  = nan(size(P,1),size(rates,2));
    for k = 1:length(params.obs_t)
        sumj(k,:) = reshape(sum(dPdrsum(k,params.r_ct{k}+1,:)./Psum(k,params.r_ct{k}+1),2),[1,size(rates,2)]);
    end

    V_mS  = -sum(sumj,1);

catch
    disp('ode15s warning ignored, sample rejected')
    V_mS = nan;
end
end