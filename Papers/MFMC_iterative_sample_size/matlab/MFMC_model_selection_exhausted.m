function [ind_for_model,xi_star] = MFMC_model_selection_exhausted(rho,C)
% function [ind_for_model,xi_star] = MFMC_model_selection(rho,C)
% MFMC model selection using exhausted method that satisfies two ordering 
% relation required by the MFMC theorem & gives the smallest objective value
%
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
% Last modified : Nov-28, 2024



% sort rho into strictly descending sequence, order saves index mapping
rho = abs(rho);
[rho_s,order] = sort(rho,'descend');

% remove same |rho|
find_zero_ind = find(rho_s(1:end-1) - rho_s(2:end)==0);
find_same_rho_ind = order(find_zero_ind+1);
order = order(~ismember(order, find_same_rho_ind));
rho = rho(order);
C = C(order);

ind = 1:length(rho);
xi_star = 1;
ind_for_model = 1;

idxs = logical(dec2bin(0:2^length(ind)-1) - '0'); % generate all 2^len(rho) subsets using 0-1 form
idxs = idxs(idxs(:,1)==1,:); % subset including the hfm with the index 1.

for i = 1:length(idxs)
    subset = ind(idxs(i,:)); % increasing index
    % disp(subset);

    if isscalar(subset) %length(subset)==1
        xi = 1;
    else
        k = subset(1:end-1);
        k_p1 = subset(2:end);
        rho_k_p2 = [rho(subset(3:end)),0];

        if  sum(C(k)./C(k_p1) <= (rho(k).^2 - rho(k_p1).^2)./(rho(k_p1).^2 - rho_k_p2.^2))==0
            kk = subset(1:end);
            rho_kk_p1 = [rho(subset(2:end)),0];
            xi = (sum(sqrt(C(kk).*(rho(kk).^2 - rho_kk_p1.^2))))^2/C(1); % add remaining elements
        else
            continue
        end
    end

    if xi <= xi_star
        ind_for_model = subset;
        xi_star = xi;
    end
end

% map the indices back
ind_for_model = order(ind_for_model);











%% BElOW are equivalent version

% function [ind_for_model,xi_star] = MFMC_model_selection_exhausted(rho,C)
% % function ind_for_model = MFMC_model_selection(rho,C)
% % select the model for MFMC that gives the smallest objective value
% % INPUT:
% % rho:      vector. correlation coefficients.
% % C:        vector. cost per sample.
% %
% % OUTPUT:
% % ind_for_model:    vector. the indices of the selected models.
% % xi_star:          scalar. value of cost efficiency.
% %
% % Nov-28, 2024
% %
%
%
% % sort rho into descending sequence, order saves index mapping
% [rho,order] = sort(rho,'descend');
% C = C(order);
%
% ind = 1:length(rho);
% xi_star = 1;
% ind_for_model = 1;
%
% idxs = logical(dec2bin(0:2^length(ind)-1) - '0'); % generate all 2^len(rho) subsets using 0-1 form
% idxs = idxs(idxs(:,1)==1,:); % subset including the hfm with the index 1.
%
%
% for i = 1:length(idxs)
%     break_flag = false;
%     subset = ind(idxs(i,:)); % increasing index
%     % disp(subset);
%
%     if length(subset)~=1
%         w = sqrt(C(subset(1))*(rho(subset(1))^2- rho(subset(2))^2)); % add first element to the objective function
%     else
%         w = sqrt(C(1));
%     end
%
%     for j=1:length(subset)-1
%         k = subset(j);
%         k_p1 =subset(j+1);
%
%         if j==length(subset)-1
%             rho_k_p2 = 0;
%         else
%             rho_k_p2 = rho(subset(j+2));
%         end
%         if   C(k)/C(k_p1) > (rho(k)^2 - rho(k_p1)^2)/(rho(k_p1)^2 - rho_k_p2^2)
%             w = w + sqrt(C(k_p1)*(rho(k_p1)^2 - rho_k_p2^2)); % add remaining elements
%         else
%             break_flag = true;
%             break
%         end
%
%     end
%     if ~break_flag
%         xi = w^2/C(1);
%         if xi<xi_star
%             ind_for_model = subset;
%             xi_star = xi;
%             disp(subset);
%             disp(xi_star);
%         end
%     end
% end
%
% % map the indices back
% ind_for_model = order(ind_for_model);