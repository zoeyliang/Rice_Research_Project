function [ind_for_model, xi_star] = MFMC_model_selection_backtrack(rho, C)

[rho, order] = sort(rho, 'descend');
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

            prev_idx = curr_indices(end);
            rho_prev = rho(prev_idx); %previous rho
            C_prev = C(prev_idx);

            rho_k = rho(k); %current rho
            C_k = C(k);

            rho_next = 0; %next rho

            % conditions to satisfy
            if (C_prev/C_k) <= (rho_prev^2 - rho_k^2)/(rho_k^2 - rho_next^2)
                continue;
            end

            % update new objective
            rho_k_vec = [rho(curr_indices),rho_k];
            rho_k_vec_next = [rho(curr_indices(2:end)),rho_k,rho_next];
            C_k_vec = [C(curr_indices),C_k];
            xi = 1/C_k_vec(1)*(sum(sqrt(C_k_vec.*(rho_k_vec.^2 - rho_k_vec_next.^2))))^2;

            % pruning
            if xi >= xi_star
                continue;
            end

            backtrack([curr_indices, k], xi, k+1);
        end
    end
end