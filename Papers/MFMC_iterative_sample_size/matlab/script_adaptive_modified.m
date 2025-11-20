% =========================================================================
% This script performs a adaptive strategy to obtain integer-valued
% sample size, with integer is selected by either FLOOR/CEIL that minimize
% the following variance:
%
%       f(x) = delta_1/m_1 + ... + delta_{k-1}/m_{k-1} + deltak/x + ...
%               sum_{j=k+1}^K(sqrt(C_j*delta_j))^2/(b_k-C_k*x)
%
%       where  b_k = b-sum_{j=1}^{k-1} C_j*m_j
%
% relative to the continuous (real-valued) optimal allocation for the MFMC
% method under a fixed computational budget p.
% -------------------------------------------------------------------------
%         |=============================|
%         |                             |
%         |   /.  .\           /@  @\   |
% Author: |     l      =-->      L      |
%         |      o                X     |
%         |   BEFORE            AFTER   |
%         |=============================|_{*Jiaxing Liang*}
%
% -------------------------------------------------------------------------
% |¬_¬|~_~|@_@|*_*|._.|-_-|G_G|×_×|/_\|w_w|Q_Q|>_>|=_=|+_+|#_#|z_z|>_<|6_9|
% -------------------------------------------------------------------------
% Last modified: Oct-22,2025
% =========================================================================

clear all
close all

%% Load data
% Example 1
% rho = [1,   9.9977e-01   9.9925e-01  9.9728e-01   9.8390e-01]; %rho_k
% rho_p1 = [rho(2:end),0];%rho_{k}
% C = [7.30e+01,7.0318e-03,1.4018e-03,5.0613e-04,2.6803e-04]; %cost
% p=828; %248 total cost


% % Example 2
% rho = [1, 9.999882e-01, 9.999743e-01, 9.958253e-01];
% rho_p1 = [rho(2:end),0];
% C=[44.395, 6.8409e-01, 2.9937e-01, 1.9908e-04];
% p=46; % total cost


% % Example 3
rho = [1, 0.9819, 0.9708];
rho_p1 = [rho(2:end),0];
C=[240.5123, 0.0166, 0.0017];
p=960;

%% ========================================================================
if p < sum(C)
    error('Budget p is insufficient: must satisfy p >= sum(C).');
end

delta = rho.^2-rho_p1.^2;
% E_real = 1/p^2*(sum(sqrt(C.*delta)))^2;

% =====================================================================
% --- Real-valued sample size: formula --- Peherstorfer
% =====================================================================
N_star = sqrt(delta./C)*p/sum(sqrt(C.*delta));
% Real-valued
f_real_value = sum(delta./N_star);
Cost_real_value = sum(C.*N_star);

% =====================================================================
% --- Initialization for possible restarts ---
% =====================================================================
restart_flag = true;

while restart_flag
    restart_flag = false;  % reset restart flag
    N_iterative = zeros(size(rho));


    for k=1:length(N_iterative)

        % =================================================================
        % --- If current samle size is floor  ---
        % =================================================================
        N_iterative_floor = floor(sqrt(delta(k)./C(k))*(p-sum(C(1:k-1).*N_iterative(1:k-1)))...
            /sum(sqrt(C(k:end).*delta(k:end)))); % remove the cost of the sample with size 1, and resample.
        N_iterative_floor = max(N_iterative_floor,1); % modify the sample size in case of 0.
        Future_re_budget_floor = p-sum(C(1:k-1).*N_iterative(1:k-1))-C(k)*N_iterative_floor;
        % E_floor = delta(k)./C(k)./N_iterative_floor/(N_iterative_floor+1);
        V_floor = sum(delta(1:k-1)./N_iterative(1:k-1)) + delta(k)/N_iterative_floor + ... %previous variance with integer + current variance with floor
            (sum(sqrt(delta(k+1:end).*C(k+1:end))))^2/Future_re_budget_floor; %Future variance with remaining budget


        % =================================================================
        % --- If current samle size is ceil  ---
        % =================================================================
        N_iterative_ceil = ceil(sqrt(delta(k)./C(k))*(p-sum(C(1:k-1).*N_iterative(1:k-1)))...
            /sum(sqrt(C(k:end).*delta(k:end)))); % remove the cost of the sample with size 1, and resample.
        Future_re_budget_ceil = p-sum(C(1:k-1).*N_iterative(1:k-1))-C(k)*N_iterative_ceil;
        % E_ceil = delta(k)./C(k)./N_iterative_ceil./(N_iterative_ceil+1);
        V_ceil = sum(delta(1:k-1)./N_iterative(1:k-1)) + delta(k)/N_iterative_ceil + ... %previous variance with integer + current variance with ceil
            (sum(sqrt(delta(k+1:end).*C(k+1:end))))^2/Future_re_budget_ceil; %Future variance with remaining budget

        % =================================================================
        % represent N_iterative(k-1) in case k=1
        % =================================================================
        if k==1
            N_iterative_k_min_1 = 0;
        else
            N_iterative_k_min_1 = N_iterative(k-1);
        end

        % =================================================================
        % Sample size non-decreasing
        % =================================================================
        N_max = max(N_iterative_ceil,N_iterative_k_min_1);

        % =================================================================
        % choose the integer with small variance
        % =================================================================
        if V_floor > V_ceil && sum(N_iterative(1:k-1).*C(1:k-1)) + sum(C(k:end))*N_max<=p
            N_iterative(k) = N_max;
        elseif N_iterative_floor < N_iterative_k_min_1 ... % this means N_iterative(k-1) >= N_iterative_ceil
                && sum(N_iterative(1:k-1).*C(1:k-1)) + sum(C(k:end))*N_iterative_k_min_1<=p
            N_iterative(k) = N_iterative_k_min_1;
        else
            N_iterative(k) = N_iterative_floor;
        end

        % =================================================================
        % --- Discard model if sample size are the same ---
        % =================================================================
        if N_iterative(k) == N_iterative_k_min_1
            delta_k_1 = rho(k-1)^2 - rho(k+1)^2;
            delta_k = rho(k+1)^2 - rho_p1(k+1)^2;

            if delta_k_1/C(k-1) < delta_k/C(k+1)
                fprintf('Discarding model index: %d\n', k);

                % Update rho, delta, and C
                rho(k)   = [];
                C(k)     = [];
                rho_p1   = [rho(2:end), 0];
                delta    = rho.^2 - rho_p1.^2;
                % Restart computation with updated variables
                restart_flag = true;
                break; % break out of k-loop to restart while-loop
            end
        end
    end
end

f_adaptive = sum(delta./N_iterative);
Cost_adaptive = sum(C.*N_iterative);

%% ========================================================================

% N_star separated by comma
formatted_elements = arrayfun(@(x) sprintf('%.2f', x), N_star, 'UniformOutput', false);
comma_separated_N_star = strjoin(formatted_elements, ', ');
formatted_elements = arrayfun(@(x) sprintf('%d', x), N_iterative, 'UniformOutput', false);
comma_separated_N_iterative = strjoin(formatted_elements, ', ');

fprintf('The real valued sample size is: [%s]\n', comma_separated_N_star);
fprintf('The real valued variance is: %f\n', f_real_value);
fprintf('The real valued cost is: %f\n', Cost_real_value);

fprintf('The adaptive sample size is: [%s]\n', comma_separated_N_iterative);
fprintf('The adaptive variance is: %f\n', f_adaptive);
fprintf('The adaptive cost is: %f\n', Cost_adaptive);










% if k <=length(N_iterative)-1
%     E_future_floor = (sum(sqrt(C(k+1:end).*delta(k+1:end))))^2/Future_re_budget_floor^2;
%     E_future_ceil = (sum(sqrt(C(k+1:end).*delta(k+1:end))))^2/Future_re_budget_ceil^2;
%
%
%     % if abs(E_future_floor-E_real)/E_future_floor > abs(E_real-E_future_ceil)/E_future_ceil && Future_re_budget_ceil>=0
%     % if E_ceil > E_floor  && Future_re_budget_ceil>=0
%     % if abs(E_floor-E_real) /abs(E_future_floor-E_real)> abs(E_ceil-E_real)/abs(E_real-E_future_ceil) && Future_re_budget_ceil>=0
%     if abs(E_future_floor-E_floor) > abs(E_ceil-E_future_ceil) && Future_re_budget_ceil>=0
%     % if abs(E_future_floor-E_real) > abs(E_real-E_future_ceil) && Future_re_budget_ceil>=0
%         N_iterative(k) = N_iterative_ceil;
%     else
%         N_iterative(k) = N_iterative_floor;
%     end
% else
%         N_iterative(k) = N_iterative_floor;
% end