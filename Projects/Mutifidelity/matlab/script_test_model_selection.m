rho = sort(rand(1,15), 'descend');
rho([5,10]) = [0.01, 0.8];
C = 10.^randi([0,3],1,15);C(9) = 0.1;

tic,[ind_for_model,xi_min] = MFMC_model_selection_exhausted(rho,C);toc
tic,[ind_for_model,xi_min] = MFMC_model_selection_backtrack(rho,C);toc


rho = sort(rand(1,10), 'descend');
C = 10 + 80*rand(1,10);