function export_chain(chain, samp_rate , nameit)

maxsample   = (length(chain.i));
minsample   = 1;
%mem         = chain.sizeGB;
params      = chain.params;

% for samplenum   = 1:fix(1/samp_rate): maxsample
%     
% rates(samplenum,:)         = chain.rates(:,samplenum);
% swap_acc(samplenum)        = chain.swap_acceptance(samplenum);
% L_post(samplenum)          = chain.MAP(samplenum,:);
% L_likel(samplenum)         = chain.MAP_LL(samplenum,:);
% omega_hist(samplenum)       = chain.omega_history(samplenum);
% 
% end 

save(['../results/',nameit,'.mat'],'chain','-v7.3');
%v7-3 is for data >2gb
disp([nameit,'.mat','--missionaccomplished']);

