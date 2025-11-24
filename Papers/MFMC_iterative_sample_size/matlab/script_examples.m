% =========================================================================
% This script compares several integer-valued MFMC sample allocation
% strategies under a fixed computational budget p, evaluated across five
% benchmark problems. The methods examined are:
%
%   • Naive floor allocation (Peherstorfer)
%   • Modified ceil–floor allocation (Gruber)
%   • Iterative allocation with modified rounding and duplicate removal (Liang)
%   • Adaptive allocation using both ceil and floor variants (Liang)
%
% -------------------------------------------------------------------------
% The five test examples are taken from the following published works:
%
% Example 1:
%   J. Liang et al.,
%   "Multilevel Monte Carlo methods for the Grad–Shafranov free boundary problem,"
%   DOI: 10.1016/j.cpc.2024.109099
%
% Example 2:
%   P. Peherstorfer et al.,
%   "Optimal model management for multifidelity Monte Carlo estimation,"
%   DOI: 10.1137/15M1046472
%
% Example 3:
%   J. Konrad et al.,
%   "Data-driven low-fidelity models for multifidelity Monte Carlo sampling
%    in plasma micro-turbulence analysis,"
%   DOI: 10.1016/j.jcp.2021.110898
%
% Example 4:
%   E. Qian et al.,
%   "Multifidelity Monte Carlo estimation of variance and sensitivity indices,"
%   DOI: 10.1137/17M1151006
%
% Example 5:
%   A. Gorodetsky et al.,
%   "A generalized approximate control variate framework for multifidelity
%    uncertainty quantification,"
%   DOI: 10.1016/j.jcp.2020.109257
%
% =========================================================================
%
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
% Last modified: Nov-23,2025
% =========================================================================

clear all
close all

%% Load data
Example = 'Example_5';

switch Example
    case 'Example_1'
        % Jiaxing's paper:
        % Multilevel Monte Carlo methods for the Grad-Shafranov free boundary problem
        % DOI: 10.1016/j.cpc.2024.109099

        rho = [1,   9.9977e-01   9.9925e-01  9.9728e-01   9.8390e-01]; %rho_k
        C = [7.30e+01,7.0318e-03,1.4018e-03,5.0613e-04,2.6803e-04]; %cost
        p=248;%248; % total cost

    case 'Example_2'
        % Peherstopher's paper:
        % Optimal model management for multifidelity Monte Carlo estimation
        % DOI: 10.1137/15M1046472

        rho = [1, 9.999882e-01, 9.999743e-01, 9.958253e-01];
        C=[44.395, 6.8409e-01, 2.9937e-01, 1.9908e-04];
        p=248;%1e4; % total cost

    case 'Example_3'
        % Konrad, Julia's paper:
        % Data-driven low-fidelity models for multi-fidelity Monte Carlo
        % sampling in plasma micro-turbulence analysis
        % DOI: 10.1016/j.jcp.2021.110898

        rho = [1, 0.9819, 0.9708];
        C=[240.5123, 0.0166, 0.0017];
        p=960;

    case 'Example_4'
        % Qian, E:
        % Multifidelity Monte Carlo estimation of variance and sensitivity indices
        % DOI: 10.1137/17M1151006

        rho=[1,0.9997,0.9465];
        C=[1, 0.05,  0.001];
        p=2000;

    case 'Example_5'
        % Gorodetsky, Alex:
        % A generalized approximate control variate framework for
        % multifidelity uncertainty quantification
        % DOI: 10.1016/j.jcp.2020.109257

        rho=[1, 0.99838, 0.99245, 0.96560, 0.70267];
        C=[1.000, 0.147, 0.026, 0.009, 0.002];
        p=955.1;
end


rho_p1 = [rho(2:end),0];
delta = rho.^2-rho_p1.^2; % dalta = rho_k^2-rho_{k+1}^2




% =====================================================================
% --- Feasibility check ---
% =====================================================================
if p < sum(C)
    error('Budget p is insufficient: must satisfy p >= sum(C).');
end



% =====================================================================
% --- Real-valued sample size: formula --- Peherstorfer
% =====================================================================
N_star = sqrt(delta./C)*p/sum(sqrt(C.*delta));
% Real-valued
f_real_value = sum(delta./N_star);
Cost_real_value = sum(C.*N_star);
% upper bound of relative difference of f
f_diff_upper_bound = sum(delta./(N_star-1))-f_real_value;

% =====================================================================
% --- Integer sample size: Naive floor --- Peherstorfer
% =====================================================================
N_floor = floor(N_star);
f_floor = sum(delta./N_floor);
Cost_floor = sum(C.*N_floor);

% =====================================================================
% --- Modified integer sample size: Ceil & floor --- Gruber
% =====================================================================
j=1;
N_floor_temp = N_floor;
N_floor_modified = N_floor_temp;
while N_floor_temp(1)<1
    N_floor_modified(j) = 1; % modify the corresponding entries
    N_floor_temp = floor(sqrt(delta(j+1:end)./C(j+1:end))*(p-sum(C(1:j).*N_floor_modified(1:j)))...
        /sum(sqrt(C(j+1:end).*delta(j+1:end)))); % remove the cost of the sample with size 1, and resample.
    j = j + 1;
end
% the remaining sample size will be calculate using the close-form solution formula.
N_floor_modified(j:end) = N_floor_temp;

f_floor_modified = sum(delta./N_floor_modified);
Cost_floor_modified = sum(C.*N_floor_modified);

% =====================================================================
% --- Integer-valued sample size: Ceil & floor --- Liang
% Iterative method with modified
% =====================================================================


% --- Initialization of possible restarts if equality of sample size occur ---
delete_equal_flag = true;
restart_flag = true;

while restart_flag
    restart_flag = false;  % reset restart flag

    N_iterative = zeros(size(rho));
    N_iterative(1) = max(floor(N_star(1)),1); % ceil if real-valued sample size falls below 1


    for k=2:length(N_iterative)

        N_iterative(k) = floor(sqrt(delta(k)./C(k))*(p-sum(C(1:k-1).*N_iterative(1:k-1)))...
            /sum(sqrt(C(k:end).*delta(k:end))));

        N_iterative(k) = max(N_iterative(k),1); % modify the sample size in case of 0.

        % --- Discard model START ---
        % Discard model when equality of sample size happens
        if delete_equal_flag
            if N_iterative(k) == N_iterative(k-1)
                delta_k_1 = rho(k-1)^2 - rho(k+1)^2; % if equal, jump to the next model
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
        % --- Discard model END ---
    end
end

f_iterative = sum(delta./N_iterative);
Cost_iterative = sum(C.*N_iterative);




%==========================================================================
% --- Integer-valued sample size: Adaptive --- Liang
%==========================================================================



% --- Initialization of possible restarts if equality of sample size occur ---
restart_flag = true;

while restart_flag
    restart_flag = false;  % reset restart flag
    N_adapt = zeros(size(rho));


    for k=1:length(N_adapt)

        % =================================================================
        % --- If current samle size is floor  ---
        % =================================================================
        N_adapt_floor = floor(sqrt(delta(k)./C(k))*(p-sum(C(1:k-1).*N_adapt(1:k-1)))...
            /sum(sqrt(C(k:end).*delta(k:end)))); % remove the cost of the sample with size 1, and resample.
        N_adapt_floor = max(N_adapt_floor,1); % modify the sample size in case of 0.
        Future_re_budget_floor = p-sum(C(1:k-1).*N_adapt(1:k-1))-C(k)*N_adapt_floor;
        % E_floor = delta(k)./C(k)./N_adapt_floor/(N_adapt_floor+1);
        V_floor = sum(delta(1:k-1)./N_adapt(1:k-1)) + delta(k)/N_adapt_floor + ... %previous variance with integer + current variance with floor
            (sum(sqrt(delta(k+1:end).*C(k+1:end))))^2/Future_re_budget_floor; %Future variance with remaining budget


        % =================================================================
        % --- If current samle size is ceil  ---
        % =================================================================
        N_adapt_ceil = ceil(sqrt(delta(k)./C(k))*(p-sum(C(1:k-1).*N_adapt(1:k-1)))...
            /sum(sqrt(C(k:end).*delta(k:end)))); % remove the cost of the sample with size 1, and resample.
        Future_re_budget_ceil = p-sum(C(1:k-1).*N_adapt(1:k-1))-C(k)*N_adapt_ceil;
        % E_ceil = delta(k)./C(k)./N_adapt_ceil./(N_adapt_ceil+1);
        V_ceil = sum(delta(1:k-1)./N_adapt(1:k-1)) + delta(k)/N_adapt_ceil + ... %previous variance with integer + current variance with ceil
            (sum(sqrt(delta(k+1:end).*C(k+1:end))))^2/Future_re_budget_ceil; %Future variance with remaining budget

        % =================================================================
        % represent N_iterative(k-1) in case k=1
        % =================================================================
        if k==1
            N_adapt_k_min_1 = 0;
        else
            N_adapt_k_min_1 = N_adapt(k-1);
        end

        % =================================================================
        % Sample size non-decreasing
        % =================================================================
        N_max = max(N_adapt_ceil,N_adapt_k_min_1);

        % =================================================================
        % choose the integer with small variance
        % =================================================================
        if V_floor > V_ceil && sum(N_adapt(1:k-1).*C(1:k-1)) + sum(C(k:end))*N_max<=p
            N_adapt(k) = N_max;
        elseif N_adapt_floor < N_adapt_k_min_1 ... % this means N_adapt(k-1) >= N_adapt_ceil
                && sum(N_adapt(1:k-1).*C(1:k-1)) + sum(C(k:end))*N_adapt_k_min_1<=p
            N_adapt(k) = N_adapt_k_min_1;
        else
            N_adapt(k) = N_adapt_floor;
        end

        % --- Discard model START ---
        % Discard model when equality of sample size happens
        if N_adapt(k) == N_adapt_k_min_1
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
        % --- Discard model END ---
    end
end

f_adaptive = sum(delta./N_adapt);
Cost_adaptive = sum(C.*N_adapt);


% =================================================================
% --- Output results ---
% =================================================================
% real-valued solution
formatted_elements = arrayfun(@(x) sprintf('%.2f', x), N_star, 'UniformOutput', false);
comma_separated_N_star = strjoin(formatted_elements, ', '); % Separated by comma
fprintf('The real valued sample size is: [%s]\n', comma_separated_N_star);
fprintf('The real valued variance is: %f\n', f_real_value);
fprintf('The real valued cost is: %f\n', Cost_real_value);
fprintf('\n');

% floor
formatted_elements = arrayfun(@(x) sprintf('%d', x), N_floor, 'UniformOutput', false);
comma_separated_N_floor = strjoin(formatted_elements, ', ');
fprintf('The floor sample size is: [%s]\n', comma_separated_N_floor);
fprintf('The floor variance is: %f\n', f_floor);
fprintf('The floor cost is: %f\n', Cost_floor);
fprintf('\n');

% modified
formatted_elements = arrayfun(@(x) sprintf('%d', x), N_floor_modified, 'UniformOutput', false);
comma_separated_N_floor_modified = strjoin(formatted_elements, ', ');
fprintf('The modified sample size is: [%s]\n', comma_separated_N_floor_modified);
fprintf('The modified variance is: %f\n', f_floor_modified);
fprintf('The modified cost is: %f\n', Cost_floor_modified);
fprintf('\n');

% iterative
formatted_elements = arrayfun(@(x) sprintf('%d', x), N_iterative, 'UniformOutput', false);
comma_separated_N_iterative = strjoin(formatted_elements, ', ');
fprintf('The iterative sample size is: [%s]\n', comma_separated_N_iterative);
fprintf('The iterative variance is: %f\n', f_iterative);
fprintf('The iterative cost is: %f\n', Cost_iterative);
fprintf('\n');

% adaptive
formatted_elements = arrayfun(@(x) sprintf('%d', x), N_adapt, 'UniformOutput', false);
comma_separated_N_adapt = strjoin(formatted_elements, ', ');
fprintf('The adaptive sample size is: [%s]\n', comma_separated_N_adapt);
fprintf('The adaptive variance is: %f\n', f_adaptive);
fprintf('The adaptive cost is: %f\n', Cost_adaptive);