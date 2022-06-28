close all; clear variables;
%%
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 14);

addpath ../sens_analysis_1.2/sens_analysis_1.2/
addpath ./functions
%% Import Data
% THIS SECTION CONTAINS ALL USER-CONTROLLED PARAMETERS
runID         = 'J50'; %indicate result ID    
nchains = 6;           %indicate number of PT chains

%%
% import chains
for i = 1:nchains                                                               
    try                                                                         
        load(['../results/',runID,'/','PTChain',num2str(i),'_cont.mat']);        
    catch                                                                       
        load(['../results/',runID,'/','PTChain',num2str(i),'.mat']);             
    end                                                                         
  chains(i) = chain;                                                            
  clear chain                                                                   
end 

chain_ct = chains;   
swap_idx = 25;


for i = 1:niterations


chain_ct = chainer_main(chain_ct,swap_idx,nchains,[],true,false );

if ~exist(['../results/',num2str(runID)],'dir')
    mkdir(['../results/',num2str(runID)])
end

for chain_idx = 1:nchains
export_chain(chain_ct(chain_idx), 1 , [char(runID),'/','PTChain',num2str(chain_idx)]);
end

end
