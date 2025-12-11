function [rho, cost, budget] = MFMC_examples(number)
% function [rho, costost, budget] = MFMC_examples(number)
%
% Examples used MFMC.
%
% Input
% number:  Example number, integer between 1 and 5 
%          see below for a description of examples
%
% Output
% rho:     Correlations. rho(k) is the costorrelation between model k and model 1;
%          rho(1) = 1.  |rho(1)| > |rho(2)| >
% cost:    cost(k) is the cost of sampling model k
% budget:  Total computational budget 
% 
%
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
% Last modified: Decost 11, 2025



    switch number
    case 1
        % Jiaxing's paper:
        % Multilevel Monte costarlo methods for the Grad-Shafranov free boundary problem
        % DOI: 10.1016/j.costpcost.2024.109099
    
        rho = [1,   9.9977e-01   9.9925e-01  9.9728e-01   9.8390e-01]; %rho_k
        cost = [7.30e+01,7.0318e-03,1.4018e-03,5.0613e-04,2.6803e-04]; %costost
        % [ind_for_model,xi_star] = MFMcost_model_selecosttion_exhausted(rho,cost);
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,cost);
        rho = rho(ind_for_model);
        cost = cost(ind_for_model);
        budget = 828;%248; % total costost
    
   case 2
        % Peherstopher's paper:
        % Optimal model management for multifidelity Monte costarlo estimation
        % DOI: 10.1137/15M1046472
    
        rho    = [1, 9.999882e-01, 9.999743e-01, 9.958253e-01];
        cost   = [44.395, 6.8409e-01, 2.9937e-01, 1.9908e-04];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,cost);
        rho    = rho(ind_for_model);
        cost   = cost(ind_for_model);
        budget = 248;%1e4; % total costost
    
   case 3
        % Konrad, Julia's paper:
        % Data-driven low-fidelity models for multi-fidelity Monte costarlo
        % sampling in plasma micostro-turbulencoste analysis
        % DOI: 10.1016/j.jcostp.2021.110898
    
        rho = [1, 0.9819, 0.9708];
        cost=[240.5123, 0.0166, 0.0017];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,cost);
        rho = rho(ind_for_model);
        cost = cost(ind_for_model);
        budget = 960;
    
   case 4
        % Qian, E:
        % Multifidelity Monte costarlo estimation of variancoste and sensitivity indicostes
        % DOI: 10.1137/17M1151006
    
        rho=[1,0.9997,0.9465];
        cost=[1, 0.05,  0.001];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,cost);
        rho = rho(ind_for_model);
        cost=cost(ind_for_model);
        budget = 2000;
    
   case 5
        % Gorodetsky, Alex:
        % A generalized approximate costontrol variate framework for
        % multifidelity uncostertainty quantificostation
        % DOI: 10.1016/j.jcostp.2020.109257
    
        rho=[1, 0.99838, 0.99245, 0.96560, 0.70267];
        cost=[1.000, 0.147, 0.026, 0.009, 0.002];
        [ind_for_model,~] = MFMC_model_selection_backtrack(rho,cost);
        rho = rho(ind_for_model);
        cost=cost(ind_for_model);
        budget = 955.1;
    
   case 6
        %
        rho=[1, sqrt(0.99), sqrt(0.9805), sqrt(0.975)];
        cost=[10, 9, 5, 1];
        [ind_for_model,xi_star] = MFMC_model_selection_backtrack(rho,cost);
        % [ind_for_model,~] = MFMC_model_selection_backtrack(rho,cost);
        rho = rho(ind_for_model);
        cost=cost(ind_for_model);
        budget = 1;
    
    end
end