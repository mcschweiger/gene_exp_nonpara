function params = chainer_init_params(opts)
% Initialize parameters object for MCMC chain from opts object

%% 
params         = opts;
params.phi     = [opts.phi_1; opts.phi_2; opts.phi_3];
params.psi     = [opts.psi_1; opts.psi_2; opts.psi_3];

params.rep_MH        = opts.rep   ;
params.rep_HMC       = opts.rep   ;
params.rep_loads     = opts.rep   ;
params.rep_q         = opts.rep   ;
params.rep_sigma     = opts.rep   ;
params.rep           = opts.rep   ;

params.i_skip        =   1        ;
