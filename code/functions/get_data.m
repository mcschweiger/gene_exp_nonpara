function [opts,runID] = get_data(dataPath,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   get_data(dataPath): fetches data from the pathway '../Datasets'
%   get_data(dataPath,Init): fetched data when initial RNA population in
%   ALL cells is known to equal Init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runID = ['1',dataPath];

vars = {'r_ct','data'};load(['../Datasets/',dataPath,'.mat'],vars{:})

%% Ensure data is a a cell, with each element a vector of RNA counts at a
%different time point
if ~iscell(r_ct)
    opts.r_ct = cell(1,size(r_ct,1));
    for i = 1:size(r_ct,1)
        opts.r_ct{i}   = r_ct(i,:);
    end
else
    opts.r_ct = r_ct;
end


%% Set parameters from data file, to be used later
% %--------------------------------------------------------------------------
if ~isfield(data,'obs_t')               %If experiment times aren't defined, assume they're linearly spaced
    data.obs_t = linspace(0,data.t_f,data.K);
end

opts.M              = max(max(r_ct));           %Assumed maximum number of RNA per cell
opts.t_0            = 0;                %First snapshot time
opts.t_f            = data.obs_t(end);  %Last snapshot time
opts.obs_t          = data.obs_t;       %List snapshots times
opts.Jk             = data.J;           %Number of cells/time
opts.K              = data.K;           %Number of times

%% Set initial probability vector
if ~isempty(varargin) %If the initial popultion is known to be zero, set Pm_init accordingly
    opts.Pm_init    = zeros(opts.M+1,1);
    opts.Pm_init(varargin{1}+1) = 1;
else                  %Otherwise, fit observed initial population to a gamma PDF and normalize
    try_hh = histcounts(opts.r_ct{1},'BinMethod','integers','Normalization','probability');
    opts.Pm_init = gampdf(0:opts.M,1,1/try_hh.Values(1))/sum(gampdf(0:opts.M,1,1/try_hh.Values(1)));
end
end

