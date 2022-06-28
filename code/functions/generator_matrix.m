function A = generator_matrix(b,rates,params)

% keyboard

rates = 10.^rates;

combos = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];
act_rates = rates(1:2*nchoosek(params.L,2)).*prod(b(combos),2);
act_rates = [act_rates(act_rates ~= 0);...
             rates(2*nchoosek(params.L,2) + find(b));...
             rates(end)];

g          = sum(b);

if g > 1
    act_combos = [nchoosek(1:g,2); flip(nchoosek(1:g,2),2)];

    for gg = 1:g    
       koff = sum(act_rates(act_combos(:,1) == gg));
       T{gg}    =  diag(act_rates(2*nchoosek(g,2)+gg)*ones(1,params.M),-1)+...
               diag(-(act_rates(2*nchoosek(g,2)+gg).*([ones(1,params.M) 0])+...
                     (0:params.M)*rates(end)+(koff)*ones(1,params.M+1)))+...
               diag((1:params.M)*rates(end),1);
    end

    A = blkdiag(T{1:g});

    for i = 1:length(act_combos(:,1))
       ii = act_combos(i,1);
       jj = act_combos(i,2);

       pos      = zeros(g,g);
       pos(ii,jj) = 1;
       kon      = diag(act_rates(prod(act_combos == [jj ii],2) ~= 0) * ones(1,size(A,1)/g));
       A        = A + kron(pos,kon);
    end
else             % Birth and Death
    A = diag(act_rates(1)*ones(1,params.M),-1)+...
        diag(-(act_rates(1).*([ones(1,params.M) 0])+(0:params.M)*act_rates(2))) +...
        diag((1:params.M)*act_rates(2),1);
end