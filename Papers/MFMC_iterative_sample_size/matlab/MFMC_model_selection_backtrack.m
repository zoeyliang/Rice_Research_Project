function [ind_for_model, xi_star] = MFMC_model_selection_backtrack(rho, C)
% function [ind_for_model,xi_star] = MFMC_model_selection(rho,C)
%
% MFMC model selection using backtrak pruning that satisfies two ordering 
% relation required by the MFMC theorem & gives the smallest objective value
% -------------------------------------------------------------------------
% INPUT
% rho           : Vector. correlation coefficients.
% C             : Vector. cost per sample.
%
% -------------------------------------------------------------------------
% OUTPUT
% ind_for_model : Vector. the indices of the selected models.
% xi_star       : Scalar. value of cost efficiency.
%
% Last modified :Nov-28, 2024



% sort rho into strictly descending sequence, order saves index mapping
rho = abs(rho);
[rho_s,order] = sort(rho,'descend');

% remove same |rho|
find_zero_ind = find(rho_s(1:end-1) - rho_s(2:end)==0);
find_same_rho_ind = order(find_zero_ind+1);
order = order(~ismember(order, find_same_rho_ind));
rho = rho(order);
C = C(order);


K = length(rho);
xi_star = 1;
global_ind = [];
current_indices = 1;  % must contain the first element

backtrack(current_indices, xi_star, 2);

ind_for_model = order(global_ind);

    function backtrack(curr_indices, curr_xi, next_k)

        % update global min
        if curr_xi <= xi_star
            xi_star = curr_xi;
            global_ind = curr_indices;
        end

        % termination
        if next_k > K
            return;
        end

        for k = next_k:K

            prev_idx = curr_indices(end); %previous index
            rho_k = rho(k); %current rho
            C_k = C(k);
            rho_next = 0; %next rho--always zero since the current is rho is the last entry

            % conditions to satisfy
            if (C(prev_idx)/C_k) <= (rho(prev_idx)^2 - rho_k^2)/(rho_k^2 - rho_next^2)
                continue;
            end

            % update new objective
            rho_k_vec = [rho(curr_indices),rho_k];
            rho_k_vec_next = [rho(curr_indices(2:end)),rho_k,rho_next];
            C_k_vec = [C(curr_indices),C_k];
            xi = 1/C(1)*(sum(sqrt(C_k_vec.*(rho_k_vec.^2 - rho_k_vec_next.^2))))^2;

            % pruning
            if xi >= xi_star
                continue;
            end

            backtrack([curr_indices, k], xi, k+1);
        end
    end
end