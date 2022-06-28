function [fixed_chain,swaps] = hungarian(chain,burn)

fixed_chain = chain;

rates     = 10.^chain.rates;

params    = chain.params;

swaps = zeros(params.L,burn-1);

indices   = [nchoosek(1:params.L,2); flip(nchoosek(1:params.L,2),2)];

for i = burn:length(chain.loglike)
    [~,anchorIdx] = max(chain.logpost);%find the index of the anchor
    anchor.loads = fixed_chain.loads(:,anchorIdx);%set up anchor
    anchor.rates = 10.^fixed_chain.rates(:,anchorIdx);
    
    temp_c.loads = chain.loads(:,i);%set up 'latest sample'
    temp_c.rates = rates(:,i);
    
    newIdx = hung_label_switching(temp_c,anchor,indices,params.L);%get switched indeces
    rateIdx = [get_t_idx(indices,params,newIdx);(newIdx+2*nchoosek(params.L,2))';size(rates,1)];
    
    fixed_chain.sigma_ast(i)    = newIdx(chain.sigma_ast(i));%switch initial gene state 
    fixed_chain.loads(:,i)      = chain.loads(newIdx,i);  %switch loads
    fixed_chain.q(:,i)          = chain.q(newIdx,i);%switch BB qs
    fixed_chain.rates(:,i)      = chain.rates(rateIdx,i);%switch production rates
    
end


function transIdx = get_t_idx(indices,params,newIdx)%This DEF works
    %obtain the correct index reordering for transition rates
    temp_indices = indices;
    transIdx     = zeros(size(indices,1),1);
    for l = 1:params.L
    temp_indices(indices == l) = newIdx(l);
    end
    lin_temp_idx = sub2ind([params.L,params.L],temp_indices(:,2),temp_indices(:,1));
    lin_idx      = sub2ind([params.L,params.L],indices(:,2),indices(:,1));
    for j = 1:size(indices,1)
    transIdx(j) = find(lin_temp_idx(j) == lin_idx);
    end
end

end