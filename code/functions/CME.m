function dP = CME(t,P,flag,rates)

global M b;

rates      = 10.^rates;
g          = sum(b);

if g > 1
    act_combos = [nchoosek(1:g,2); flip(nchoosek(1:g,2),2)];

    for gg = 1:g    
       koff = sum(rates(act_combos(:,1) == gg));
       T{gg}    =  diag(rates(2*nchoosek(g,2)+gg)*ones(1,M),-1)+...
               diag(-(rates(2*nchoosek(g,2)+gg).*([ones(1,M) 0])+...
                     (0:M)*rates(end)+(koff)*ones(1,M+1)))+...
               diag((1:M)*rates(end),1);
    end

    A = blkdiag(T{1:g});

    for i = 1:length(act_combos(:,1))
       ii = act_combos(i,1);
       jj = act_combos(i,2);

       pos      = zeros(g,g);
       pos(ii,jj) = 1;
       kon      = diag(rates(prod(act_combos == [jj ii],2) ~= 0) * ones(1,size(A,1)/g));
       A        = A + kron(pos,kon);
    end
else             % Birth and Death
    A = diag(rates(1)*ones(1,M),-1)+...
        diag(-(rates(1).*([ones(1,M) 0])+(0:M)*rates(2))) +...
        diag((1:M)*rates(2),1);
end
% keyboard

 dP = A*P;
end
