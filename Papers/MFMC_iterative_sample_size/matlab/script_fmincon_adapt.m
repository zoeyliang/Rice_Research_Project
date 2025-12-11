%% Load data
Example = 'Example_5';

 MFMC_examples(number)
 

switch Example
    case 'Example_1'
        % Jiaxing's paper:
        % Multilevel Monte Carlo methods for the Grad-Shafranov free boundary problem
        % DOI: 10.1016/j.cpc.2024.109099

        rho = [1,   9.9977e-01   9.9925e-01  9.9728e-01   9.8390e-01]; %rho_k
        C = [7.30e+01,7.0318e-03,1.4018e-03,5.0613e-04,2.6803e-04]; %cost
        % [ind_for_model,xi_star] = MFMC_model_selection_exhausted(rho,C);
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,C);
        rho = rho(ind_for_model);
        C=C(ind_for_model);
        b=828;%248; % total cost

    case 'Example_2'
        % Peherstopher's paper:
        % Optimal model management for multifidelity Monte Carlo estimation
        % DOI: 10.1137/15M1046472

        rho = [1, 9.999882e-01, 9.999743e-01, 9.958253e-01];
        C=[44.395, 6.8409e-01, 2.9937e-01, 1.9908e-04];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,C);
        rho = rho(ind_for_model);
        C=C(ind_for_model);
        b=248;%1e4; % total cost

    case 'Example_3'
        % Konrad, Julia's paper:
        % Data-driven low-fidelity models for multi-fidelity Monte Carlo
        % sampling in plasma micro-turbulence analysis
        % DOI: 10.1016/j.jcp.2021.110898

        rho = [1, 0.9819, 0.9708];
        C=[240.5123, 0.0166, 0.0017];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,C);
        rho = rho(ind_for_model);
        C=C(ind_for_model);
        b=960;

    case 'Example_4'
        % Qian, E:
        % Multifidelity Monte Carlo estimation of variance and sensitivity indices
        % DOI: 10.1137/17M1151006

        rho=[1,0.9997,0.9465];
        C=[1, 0.05,  0.001];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,C);
        rho = rho(ind_for_model);
        C=C(ind_for_model);
        b=2000;

    case 'Example_5'
        % Gorodetsky, Alex:
        % A generalized approximate control variate framework for
        % multifidelity uncertainty quantification
        % DOI: 10.1016/j.jcp.2020.109257

        rho=[1, 0.99838, 0.99245, 0.96560, 0.70267];
        C=[1.000, 0.147, 0.026, 0.009, 0.002];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,C);
        rho = rho(ind_for_model);
        C=C(ind_for_model);
        b=955.1;

      case 'Example_6'

        rho=[1, sqrt(0.99), sqrt(0.9805), sqrt(0.975)];
        C=[10, 9, 5, 1];
        [ind_for_model,xi_star] = MFMC_model_selection_backtrack(rho,C);
        % [ind_for_model,~] = MFMC_model_selection_backtrack(rho,C);
        rho = rho(ind_for_model);
        C=C(ind_for_model);
        b=1;

end
% Compute Delta_k = sigma_1^2 * (rho_{1,k}^2 - rho_{1,k+1}^2)
rho_p1 = [rho(2:end), 0];
Delta = rho.^2 - rho_p1.^2;   % Delta_k

K = length(C);

% Precompute S_j and B_j for all stages
S = zeros(1, K+1);  % S_j = sum_{k=j}^K sqrt(C_k * Delta_k), S_{K+1} = 0
B = zeros(1, K);    % B_j = sqrt(Delta_j / C_j)

% Compute S_j (backwards from K to 1)
S(K+1) = 0;
for j = K:-1:1
    S(j) = sqrt(C(j) * Delta(j)) + S(j+1);
end

% Compute B_j
for j = 1:K
    B(j) = sqrt(Delta(j) / C(j));
end

% Display precomputed values
fprintf('Precomputed values:\n');
for j = 1:K
    fprintf('S_%d = %.6f, B_%d = %.6f\n', j, S(j), j, B(j));
end
fprintf('\n');

% Preallocate
m_opt = zeros(1, K);  % Continuous optimal solution
m_int = ones(1, K);   % Integer solution
f_val = ones(1, K);   % Integer solution


% Loop over fidelities
for j = 1:K
    fprintf('\n--- Solving optimization for j = %d ---\n', j);
    [m_opt, m_int,f_val] = solve_opt(j, C, Delta, b, K, m_int, m_opt,f_val, S, B);
end

% Display final results
fprintf('\n========== Final Results ==========\n');
fprintf('Continuous optimal solution:\n');
disp(m_opt);
fprintf('Integer solution:\n');
disp(m_int);
fprintf('Objective:\n');
disp(f_val);
fprintf('Total objective value: %.6f\n', sum(Delta./m_int));
fprintf('Total budget used: %.6f\n', sum(C.*m_int));
fprintf('Remaining budget: %.6f\n', b - sum(C.*m_int));

% Verify monotonicity constraints
fprintf('\nMonotonicity check:\n');
all_monotonic = true;
for k = 2:K
    if m_int(k) < m_int(k-1)
        fprintf('  Violation at k=%d: m(%d)=%d < m(%d)=%d\n',...
            k, k, m_int(k), k-1, m_int(k-1));
        all_monotonic = false;
    end
end
if all_monotonic
    fprintf('  All monotonicity constraints satisfied.\n');
end

function [m_opt, m_int,f_val] = solve_opt(j, C, Delta, b, K, m_int, m_opt,f_val, S, B)
% Compute remaining budget b_j = b - sum_{k<j} C_k * overline_m_k
if j > 1
    b_j = b - sum(C(1:j-1) .* m_int(1:j-1));
else
    b_j = b;
end

fprintf('Remaining budget b_%d = %.6f\n', j, b_j);

% If no budget left, use previous stage value
if b_j <= 0
    fprintf('Warning: No budget remaining at stage %d\n', j);
    m_opt(j) = m_int(j-1);
    m_int(j) = m_int(j-1);
    return;
end

% Constant term from fixed previous stages
 const_term = sum(Delta(1:j-1) ./ m_int(1:j-1));



% Use precomputed S_j and B_j
S_j_plus_1 = S(j+1);  % S_{j+1} = sum_{k=j+1}^K sqrt(C_k * Delta_k)
% B_{j+1} for monotonicity constraint
if j < K
    B_j_plus_1 = B(j+1);
else
    B_j_plus_1 = 0;
end


% Define objective function F(m_j)
fobj = @(x) const_term + Delta(j)/x + (j < K) * (S_j_plus_1)^2 / (b_j - C(j)*x);

% Compute lower and upper bounds
% Lower bound: monotonicity m_j >= m_{j-1} (or 1 for first stage)
if j == 1
    lb = 1;
else
    lb = m_int(j-1);
end

% Upper bound 1: budget constraint b_j - C_j*x > 0
ub1 = b_j / C(j);

% Upper bound 2: from monotonicity constraint m_{j+1} >= m_j
% Derived from: m_{j+1} = (B_{j+1} * (b_j - C_j*m_j)) / S_{j+1} >= m_j
if j < K
    ub2 = (B_j_plus_1 * b_j) / (S_j_plus_1 + C(j) * B_j_plus_1);
    % Effective upper bound is the minimum of the two
    ub = min(ub1, ub2);

    fprintf('Upper bounds: budget constraint=%.6f, monotonicity=%.6f\n',...
        ub1, ub2);
    fprintf('S_%d = %.6f, B_%d = %.6f\n', j+1, S_j_plus_1, j+1, B_j_plus_1);
else
    ub = ub1;
end

% Check feasibility
if lb > ub
    fprintf('Warning: Infeasible bounds at stage %d (lb=%.6f > ub=%.6f)\n',...
        j, lb, ub);
    % Use upper bound as solution
    m_opt(j) = ub;
    m_int(j) = ceil(ub);  % Round up to stay within budget
    return;
end

fprintf('Bounds: lb=%.6f, ub=%.6f\n', lb, ub);

% Initial point: midpoint of bounds, ensuring feasibility
x0 = max(lb + 0.01, min((lb + ub)/2, ub - 0.01));

% Use fminbnd for single-variable optimization
options = optimset('Display', 'iter', 'TolX', 1e-12);
[x_star, fval, exitflag] = fminbnd(fobj, lb, ub, options);

fprintf('Continuous optimal m_%d = %.6f, objective = %.6f\n',...
    j, x_star, fval);

f_val(j) = fval;
% Store continuous optimal solution
m_opt(j) = x_star;

% Integer rounding with look-ahead consideration
candidates = [floor(x_star), ceil(x_star)];
best_val = Inf;
best_m = candidates(1);

for candidate = candidates
    % Check feasibility: monotonicity and budget constraints
    is_feasible = (candidate >= lb);
    
    % For all stages, budget constraint is C(j)*candidate < b_j
    % (strict inequality to ensure positive remaining budget for j < K)
    if j < K
        is_feasible = is_feasible && (C(j)*candidate < b_j);
    else
        % For j = K, use <= to allow using all remaining budget
        is_feasible = is_feasible && (C(j)*candidate <= b_j);
    end
    
    % If feasible, evaluate objective
    if is_feasible
        candidate_val = fobj(candidate);
        if candidate_val < best_val
            best_val = candidate_val;
            best_m = candidate;
        end
    end
end

% If no feasible candidate found (should not happen with proper bounds)
if isinf(best_val)
    % Fallback: use the continuous optimum or previous stage value
    if j > 1
        best_m = m_int(j-1);
    else
        best_m = max(1, floor(x_star));
    end
    fprintf('Warning: No feasible integer candidate found at stage %d, using %d\n', j, best_m);
end

m_int(j) = best_m;

fprintf('Integer m_%d = %d\n', j, m_int(j));
end







%% ========================================================================
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