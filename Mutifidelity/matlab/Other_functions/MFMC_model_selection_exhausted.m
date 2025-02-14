function [ind_for_model,sum_sqrt_C_rho] = MFMC_model_selection_exhausted(rho,C)
% function ind_for_model = MFMC_model_selection(rho,C)
% select the model for MFMC that gives the smallest objective value
% INPUT:
% rho:      vector. correlation coefficients.
% C:        vector. cost per sample.
%
% OUTPUT:
% ind_for_model:    vector. the indices of the selected models.
% sum_sqrt_C_rho:   scalar. value of \sum_k sqrt{C_k(rho_{1,k}^2-rho_{1,k+1}^2)}
%                   W^{MF} = sum_sqrt_C_rho^2.
%
% Nov-28, 2024
%


% sort rho into descending sequence, order saves index mapping
[rho,order] = sort(rho,'descend');
C = C(order);

ind = 1:length(rho);
W_star = C(1);
ind_for_model = 1;

idxs = logical(dec2bin(0:2^length(ind)-1) - '0'); % generate all 2^len(rho) subsets using 0-1 form
idxs = idxs(idxs(:,1)==1,:); % subset including the hfm with the index 1.


for i = 1:length(idxs)
    subset = ind(idxs(i,:)); % increasing index 
    % disp(subset);
    
    if length(subset)~=1
        w = sqrt(C(subset(1))*(rho(subset(1))^2- rho(subset(2))^2)); % add first element to the objective function
    else
        w = sqrt(C(1));
    end

    for j=1:length(subset)-1
        k = subset(j);
        k_p1 =subset(j+1);

        if j==length(subset)-1
            rho_k_p2 = 0;
        else
            rho_k_p2 = rho(subset(j+2));
        end
        if   C(k)/C(k_p1) > (rho(k)^2 - rho(k_p1)^2)/(rho(k_p1)^2 - rho_k_p2^2)
            w = w + sqrt(C(k_p1)*(rho(k_p1)^2 - rho_k_p2^2)); % add remaining elements
        else
            continue
        end

    end

    if w^2<W_star
        ind_for_model = subset;
        W_star = w^2;
    end
end
sum_sqrt_C_rho = sqrt(W_star);

% map the indices back
ind_for_model = order(ind_for_model);

