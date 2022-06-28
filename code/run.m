clear variables;
close all;
clc
format compact
addpath ../sens_analysis_1.2/sens_analysis_1.2
addpath ./functions

%% Import Data
% THIS SECTION CONTAINS ALL USER-CONTROLLED PARAMETERS


path         = 'J50';           %specify the name of a .mat file in the folder: '../Datasets'
[opts,runID] = get_data(path,0);%fetch data
opts.nchains            = 6;    %indicate the number of parallel chains
niterations             = 1e4;  %indicate the maximum number of MCMC iterations

%% Set parameters
% LEAVE AS DEFAULT, THESE PARAMETERS DO NOT NEED TO BE CHANGED


% Sampler parameters
opts.tag = 1; % Debugging param. leave at default

opts.L          = 8;                               % Weak limit on possible gene states

opts.n_params   = 2*nchoosek(opts.L,2) + opts.L+1; % # of Rate parameters

opts.q_joint    = 1;                               % # of success probs. per iter

opts.n_joint    = 2;                               % # rates jointly sampled together

indices         = ...                              % get names for all kinetic params.
    [nchoosek(1:opts.L,2); flip(nchoosek(1:opts.L,2),2)];
for i = 1:size(indices,1)
   knames{i}  = ['k_{' num2str(indices(i,:),'%d') '}'];
end
for i = 1:opts.L
    bnames{i} = ['beta_' num2str(i)];
end
opts.names   = [knames bnames {'gamma'}];

opts.rep       = 2; % # samples per iteration

% MH proposal dist. parameters
opts.p_lambda  =  repmat(.0001,opts.n_params,1);

opts.kappa     =  .1;

opts.cov       = opts.p_lambda.*eye(size(opts.names,2));

% Statistical parameters
% Prior parameters
opts.phi_1 = repmat(-1,size(knames,2),1);opts.phi_2 = repmat(-1,size(bnames,2),1);opts.phi_3 = -1;
opts.psi_1 = repmat(.5,size(knames,2),1);opts.psi_2 = repmat(.5,size(bnames,2),1);opts.psi_3 = .5;

% Hyperprior parameters
opts.alpha  = 1;opts.zeta   = 5;opts.xi     = .5;


%% HMC parameters
opts.HMC_L      = 20;opts.HMC_eps    = 5*10^-5;

opts.phi_m      = 8 ;opts.psi_m      = 1;

opts.relstep_init   = opts.p_lambda;


%% Parallel tempering params.
%Omega
opts.max_omega                = 1.0;

opts.omega_init               = 0.666;

opts.adapt_omega_rate         = 0.04;

opts.adapt_omega_interval     = 10;

%Step size
opts.optimal_swap_acceptance = 0.14;

opts.adapt_relstep_rate      = 0.2;

opts.adapt_relstep_interval  = 1;

opts.optimal_step_acceptance = 0.24;

opts.adapt_last              = 150; 

opts.max_init_steps          = 50;

opts.min_adaptation_factor   = 0.8;

opts.max_adaptation_factor   = 1.25;

%Settings

opts.ParTemp            = 1;

opts.param_scale        = [opts.p_lambda];

opts.nparams            = size(opts.param_scale,1);

opts.shuffle            = 1;

opts.nswaps             = 100; %10000

opts.chain_MCMCsteps    = 30;

workers                 = 6;


% Set random number seed

if (opts.shuffle)
    globalstream = RandStream('mrg32k3a','Seed','shuffle');
else
    globalstream = RandStream('mrg32k3a','Seed',1,'NumStreams');
end

RandStream.setGlobalStream( globalstream );

opts.update_stepsize_fcn = @(relstep) repmat(relstep, [1 opts.nparams]).*(opts.param_scale);

%% INITIALIZE/DISPLAY PT CHAINS

fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Initializing chains...\n');

chain_ct = chainer_main([],0,opts.nchains,opts);

fprintf(1,'------------------------------------------------------\n');

swap_idx = opts.nswaps;nchains  = opts.nchains;
clear opts

%%
for i = 1:niterations
chain_ct = chainer_main(chain_ct,swap_idx,nchains,[]);

if ~exist(['../results/',num2str(runID)],'dir')
    mkdir(['../results/',num2str(runID)])
end

for chain_idx = 1:nchains
export_chain(chain_ct(chain_idx), 1 , [char(runID),'/','PTChain',num2str(chain_idx)]);
end

end




