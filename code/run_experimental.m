clear variables;
close all;
clc
format compact

addpath ../sens_analysis_1.2/sens_analysis_1.2
addpath ./functions
%% Import Data
% THIS SECTION CONTAINS ALL USER-CONTROLLED PARAMETERS

load('../Datasets/data-fig2g.mat')
data = data_tot;         % load total numbers
opts.type   = 'f';
runID       = ['3fig2',opts.type];
niterations             = 1e4; %indicate the maximum number of iterations desired



%%
r_ct = cell(size(data));
Jk   = zeros(1,length(ts)-3);
for i = 1:size(r_ct,2)-3      % Remove NaNs
    r_ct{i}      = round(data{i});
    ind          = find(isnan(r_ct{i}));
    r_ct{i}(ind) = [];
    if ~isempty(r_ct{i})     % Find the maxmimum RNA count per timepoint
        max_rna(i)       = max(r_ct{i});
    end
    Jk(i)        = size(r_ct{i},1);
end
r_ct = r_ct(~cellfun('isempty',r_ct));

timepoints = ts(1:end-3);
if size(r_ct,2) ~= size(timepoints,2)
    disp('Warning: timepoints and samples do not match')
end

% Initial Data
f = figure('visible','off');
hh = histogram(r_ct{1},'Normalization','pdf');
hold on
x  = linspace(0,max(r_ct{1}),100);
plot(x,gampdf(x,1,1/hh.Values(1)))

%% Data Params
% %--------------------------------------------------------------------------
opts.r_ct           = r_ct;
opts.M              = max(max_rna);
opts.t_0            = timepoints(1);
opts.t_f            = timepoints(end);
opts.obs_t          = timepoints;
opts.Jk             = Jk;
opts.K              = size(r_ct,2);

opts.Pm_init = gampdf(0:opts.M,1,1/hh.Values(1))'./sum(gampdf(x,1,1/hh.Values(1)));
% opts.Pm_init    = zeros(opts.M+1,1);
% opts.Pm_init(1) = 1;


%% init chain
opts.tag = 1;

opts.L          = 6;                               % # of possible gene states
opts.n_params   = 2*nchoosek(opts.L,2) + opts.L+1; % # of Parameters
opts.n_k        = 2*nchoosek(opts.L,2);
opts.n_beta     = opts.L;
opts.L_prime    = min(opts.L,6);                   % # loads sampled per iter

opts.q_joint    = 1;

% # rates jointly sampled together
opts.n_joint    = 2;

indices   = [nchoosek(1:opts.L,2); flip(nchoosek(1:opts.L,2),2)];
for i = 1:size(indices,1)
   knames{i}  = ['k_{' num2str(indices(i,:),'%d') '}'];
end
for i = 1:opts.L
    bnames{i} = ['beta_' num2str(i)];
end
opts.names   = [knames bnames {'gamma'}];

% Prior parameters
opts.phi_1   = repmat(-1,size(knames,2),1);
opts.phi_2   = repmat(-1,size(bnames,2),1);
opts.phi_3   = -1;

opts.psi_1   = repmat(.5,size(knames,2),1);
opts.psi_2   = repmat(.5,size(bnames,2),1);
opts.psi_3   = .5;

opts.alpha  = 1;
opts.zeta   = 5;
opts.xi     = .5;

% Shape parameter for proposal
opts.p_lambda  =  repmat(.0001,opts.n_params,1);
opts.kappa     =  .1;
opts.cov       = opts.p_lambda.*eye(size(opts.names,2));

opts.rep = 2; % # samples per iteration

opts.HMC_L      = 20;
opts.HMC_eps    = 5*10^-5;
opts.phi_m      = 8;
opts.psi_m      = 1;

opts.relstep_init   = opts.p_lambda;

%% ---------DEBUGGING RATES and PROBS-----------

%% ----------PARELLEL TEMPERING PARAMETERS *OMEGA* -------

opts.max_omega                = 1.0;

opts.omega_init               = 0.666;

opts.adapt_omega_rate         = 0.04;

opts.adapt_omega_interval     = 10;

%% ----------PARELLEL TEMPERING PARAMETERS *STEP SIZE* -------

opts.optimal_swap_acceptance = 0.14;

opts.adapt_relstep_rate      = 0.2;

opts.adapt_relstep_interval  = 1;

opts.optimal_step_acceptance = 0.24;

opts.adapt_last              = 150; %2900

opts.max_init_steps          = 50; %500

opts.min_adaptation_factor   = 0.8;

opts.max_adaptation_factor   = 1.25;

%% ----------PARELLEL TEMPERING SETTINGS PARAMETERS -------

opts.ParTemp            = 1;

opts.param_scale        = [opts.p_lambda];

opts.nparams            = size(opts.param_scale,1);

opts.nchains            = 10;%3;

opts.shuffle            = 1;

opts.nswaps             = 50; %10000

opts.chain_MCMCsteps    = 30;

workers                 = 6;


% Nlabs

if (opts.shuffle)
    globalstream = RandStream('mrg32k3a','Seed','shuffle');
else
    globalstream = RandStream('mrg32k3a','Seed',1,'NumStreams');
end

RandStream.setGlobalStream( globalstream );


% update stepsize
%  template: [stepsize] = @(relstep)
opts.update_stepsize_fcn = @(relstep) repmat(relstep, [1 opts.nparams]).*(opts.param_scale);

%% ----------INITIALIZE/DISPLAY MULTIPLE CHAINS----------

fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Initializing chains...\n');



chain_ct = chainer_main([],0,opts.nchains,opts);


fprintf(1,'------------------------------------------------------\n');
% sampler_update_display( swap_idx, opts, beta, relative_step_size, step_acceptance, swap_acceptance, energy_chain )

swap_idx = opts.nswaps;
nchains  = opts.nchains;

clear opts


for i = 1:niterations


chain_ct = chainer_main(chain_ct,swap_idx,nchains);

if ~exist(['./results/',num2str(runID)],'dir')
    mkdir(['./results/',num2str(runID)])
end

for chain_idx = 1:nchains
export_chain(chain_ct(chain_idx), 1 , [runID,'/','PTChain',num2str(chain_idx)]);
end

end




