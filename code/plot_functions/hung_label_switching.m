function assignment = hung_label_switching(chain_new,chain_anchor,Idx,L)
%%% rates_new (beta_{m}^{i}, k_{m,n}^{i}) : sampled rates at
%%% the i-th iteration
%%% rates_MAP (beta_{m}^{*}, k_{m,n}^{*})  : this is the reference rates acquired based on
%%% MAP until the ith iteration
%%% cost(m,mp) = d_{m,mp} = \sum_{n=1}^{M}|(beta_m^{i})*(k_{m,n}^{i})-(beta_{mp}^{*})*(k_{mp,n}^{*})|

for i = 1:L
    outIdx(i,:) = find(Idx(:,1)==i)';%outIdx(i,:) will be indices in rates of transitions out of i
end



for m = 1:L
    beta_m =  chain_new.rates(end-L+m-1) * chain_new.loads(m) ;
    for mp = 1:L
        beta_mp =  chain_anchor.rates(end-L+mp-1) * chain_anchor.loads(mp);
        switch sum(chain_anchor.loads(:))
            case 1
                cost_matrix(m,mp) = abs(beta_m-beta_mp);
            case 2
                cost_matrix(m,mp) = abs(beta_m-beta_mp) + ...
                    10*(1-(chain_new.loads(m) == chain_anchor.loads(mp)));
            case 3
                
                cost_matrix(m,mp) = abs(beta_m-beta_mp) + ...
                    10*(1-(chain_new.loads(m) == chain_anchor.loads(mp)));
                
            case 4
                cost_matrix(m,mp) = abs(beta_m-beta_mp).^2 + ...
                    10*(1-(chain_new.loads(m) == chain_anchor.loads(mp)));          
            case 5
                cost_matrix(m,mp) = sum(abs(beta_m .* chain_new.rates(outIdx(m,:))-...
                    beta_mp .* chain_anchor.rates(outIdx(mp,:))));%kdifference
        end
        
        %         cost_matrix(m,mp) = abs(beta_m-beta_mp) + ...
        %             10*(1-(chain_new.loads(m) == chain_anchor.loads(mp)));
        
        
    end
end
% Applying the Hungarian algorithm
[assignment,~] = munkres(cost_matrix');%returns the correct permutation of

end

%         switch cost
%             case 1
%                 cost_matrix(m,mp) = abs(beta_m-beta_mp) + ...(1-(chain_new.loads(m) == chain_anchor.loads(mp)))*1e10;%difference
%             case 2
%                 cost_matrix(m,mp) = abs(beta_m-beta_mp)^2;%l2 norm
%             case 3
%                 intf = @(t) abs(gampdf(t,beta_m)-gampdf(t,beta_mp));
%                 cost_matrix(m,mp) = int(intf(x),x,0,Inf);
%             case 3

%             case 4
%                 cost_matrix(mp,m) = sum(abs(beta_m .* chain_new.rates(outIdx(m,:))-...
%                     beta_mp .* chain_anchor.rates(outIdx(mp,:))));%l2 kdifference
%         end