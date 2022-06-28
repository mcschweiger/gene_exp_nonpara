function [r_ct,data,ground,units] = generate_simulated_data(rates,...
                                                        K,...
                                                        J,...
                                                        init,...
                                                        n_sim,...
                                                        t_0,t_f,obs_t,...
                                                        G,show_data)
%% Units

t_units   = 's';            %time units [sec]
obs_units = 'cts';          %obs  units [counts]
r_units   = '1/s';          %rate units [1/sec]
r_d_units = '1/(mol*s)';    %rate degradation units [1/( sec * mol )]

%%%%%%%%% Gillepie Simulation   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rates
K_lj    = rates(1:nchoosek(G,2));                       % Low to High Gene State Rates
K_jl    = rates(nchoosek(G,2)+1:nchoosek(G,2)*2);       % High to Low Gene State Rates
Beta_l  = rates(2*nchoosek(G,2)+1:2*nchoosek(G,2)+G);   % Production Rates
gam     = rates(end);                                   % Degradation Rate

%% Updates:
V = [];
for i = 1:G-1
    V = [V [repelem(0,i-1,G-i); repelem(-1,G-i); diag(1*ones(1,G-i))]]; % Low to High Gene State Transitions
end
V     = [V -V];  % Include High to Low Gene State Transitions

V     = [[V; zeros([1 size(V,2)])] [zeros(G); ones([1 G])]]; % Production from each Gene State
V     = [V [zeros([G 1]); -1]];

%%

% obs_t = linspace(t_0,t_f,K);

events  = zeros(2*nchoosek(G,2)+G+1,n_sim);

% Initialize population arrays
pop_s        = zeros(n_sim,K,G+1);

for iter = 1:n_sim %(>=JK)
% cell trajectories
    % Initilaize the population for each iteration

    % Initialize index
    i       = 2;
    k       = obs_t(1);
    
    
    rc_init = init.rc_init;
    rc_init(end) = Discrete_sampler(init.m_init')-1;
% keyboard
    
    updates = [rc_init' obs_t(1)];
    population   = rc_init;        %column vector

        while k < obs_t(end)                    % update between timepoints
            % current population
            g = population(1:G)';    % Gene state
            m = population(G+1);     % # mRNA

            % Propensities
            g_rep_lj = [];
            g_rep_jl = [];
            for h = 1:G-1
               g_rep_lj = [g_rep_lj repelem(g(h),G-h)];
               g_rep_jl = [g_rep_jl g(h+1:G)];
            end

% keyboard
            props(1:nchoosek(G,2))                       = K_lj.*g_rep_lj; %k_lj
            props(nchoosek(G,2)+1:2*nchoosek(G,2))       = K_jl.*g_rep_jl; %k_jl
            props(2*nchoosek(G,2)+1:2*nchoosek(G,2)+G)   = Beta_l.*g;      % Produce
            props(2*nchoosek(G,2)+G+1)                   = gam*m;          % Degrade

            % Draw the event and timestep
                    q = rand;                            % Independent uniformly ..(0,1)
            prop_sum  = sum(props);                      % Sum of propensities
            rxn_probs = props./prop_sum;                 % Compute discrete probabilities for each reaction

            % Determine reaction
            j     = 1;
            p_sum = 0;

            while p_sum < q

                p_sum = p_sum+rxn_probs(j);
                j     = j+1;

            end

            reaction         = j-1;

            events(j-1,iter) = events(j-1,iter) + 1;

            % Determine timestep
            dt = exprnd(1/prop_sum);

            % Update population and time
            population = population + V(:,reaction);
            k          = k + dt;
            updates(i,:) = [population' k];
            i = i+1;
        end

        for k = 1:length(obs_t)
            above = find(updates(:,end) > obs_t(k),1);
            pop_s(iter,k,:) = updates(above - 1,1:end-1);
        end

    if mod(iter,1000) == 0
        disp(['Gillespie Iteration: ',num2str(iter)])
    end
end
%% Visualize single Trajectory
if show_data
    gene_traj = sum(updates(:,1:G) .* [1:G],2);

    figure
    subplot(1,2,1)
    stairs(updates(:,end),gene_traj,'-o')
    yticks([1:G])
    grid on
    xlabel('t (s)')
    ylabel('Gene State')
    xlim([0 updates(end,end)])
    ylim([1-.1 G+.1])
    title('Sample Gene State Trajectory')

    subplot(1,2,2)
    plot(updates(:,end),updates(:,G+1))
    xlabel('t (s)')
    ylabel('mRNA (count)')
    xlim([0 updates(end,end)])
    title('Sample mRNA State Trajectory')

    drawnow
end
%% Average # of Events

means = nan(size(events,1),1);

for i = 1:size(events,1)

    means(i) = mean(events(i,:));

end

disp( repelem('.',size(['Low G to High G: ',num2str(means(1:nchoosek(G,2))')],2)))
disp('Average Number of Events')
disp(repelem('-',size('Average Number of Events',2)))
disp(['Low G to High G: ',num2str(means(1:nchoosek(G,2))')])
disp(['High G to Low G: ',num2str(means(nchoosek(G,2)+1:2*nchoosek(G,2))')])
disp(['Productions    : ',num2str(means(2*nchoosek(G,2)+1:2*nchoosek(G,2)+G)')])
disp(['Degradations   : ',num2str(means(end))])
disp( repelem('.',size(['Low G to High G: ',num2str(means(1:nchoosek(G,2))')],2)))
%% Construct Data set
c_pool = pop_s(:,:,G+1);                             % Pull mRNA values

% Initialize Dataset

r_ct = zeros(size(obs_t,2),J);

% Sample cells for each timepoint

for k = 1:size(obs_t,2)

    if size(c_pool,1) >= J

        c_s             = randsample(size(c_pool,1),J); % Random sample of cells

        r_ct(k,:,:)     = c_pool(c_s,k,:);              % Save mRNA # for selected cells

        c_pool(c_s,:,:) = [];                           % Remove the selected cells from the pool

    else

        disp(['Maximum Number of cells reached before t =',...
                                                          num2str(obs_t(k))])
        obs_t  = obs_t(1:k);

        break

    end
end

max_m = max(max(r_ct));

disp(repelem('.',size(['Maximum mRNA: ',num2str(max_m)],2)))
disp(['Maximum mRNA: ',num2str(max_m)])
disp(repelem('.',size(['Maximum mRNA: ',num2str(max_m)],2)))
%% Ground truth

ground.rates   = rates;
ground.rc_init = rc_init;
ground.G       = G;

%% Data outputs

data.m_ev   =  means;
data.M      = max_m;
data.t_f    = t_f;
data.J      = J;
data.K      = K;
data.obs_t  = obs_t;

%% Units

units.t_units   = t_units  ;          %time units [sec]
units.obs_units = obs_units;          %obs  units [counts]
units.r_units   = r_units;          %rate units [1/sec]
units.r_d_units = r_d_units;    %rate degradation units [1/( sec * mol )]


end