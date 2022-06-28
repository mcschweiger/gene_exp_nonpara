function export_MAP(MAP,params,nameit)


maxsample   = (length(MAP));
minsample   = 1;


cd '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project2_old_HMM/Project2_HMM_part/HMM_MJP/results' ;
save([nameit,'.mat'],'params','MAP','-v7.3');
%v7-3 is for data >2gb
disp([nameit,'.mat','--missionaccomplished']);
cd '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project2_old_HMM/Project2_HMM_part/HMM_MJP' ;
