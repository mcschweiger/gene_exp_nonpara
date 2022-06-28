function [LLhood,LPost] = calculate_Lpost_est(sample,params)

%% Part1 - Contribution from the likelihood 
LLhood = -sample.l_hood;

%% Part2 - contribution from the prior
LLPriors = params.phi.*log(params.phi./params.psi) - gammaln(params.phi)...
         + (params.phi-1).*log(sample.rates) - (params.phi./params.psi).*sample.rates;

%% FINAL MAP
LPost      = LLhood + sum(-LLPriors);


