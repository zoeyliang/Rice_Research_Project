function plot_n_print(Result_test_covar)
%function plot_n_print(Result_test_covar): tables n plots
%INPUT:
%Result_test_covar: structure.
%OUTPUT:
%plot of tables etc for the input.

T = table;
T.iter = {'j-1';'j'};
T.MF_variance = [Result_test_covar.V_MF_j_1; Result_test_covar.V_MF_j];
T.efficiency_xi = [Result_test_covar.xi_j_1; Result_test_covar.xi_j];
T.sample_size = int16([Result_test_covar.N_j_1; Result_test_covar.N_j]);
T.coefficient_alpha = [Result_test_covar.alpha_j_1; Result_test_covar.alpha_j];
disp(T)

S = table;
S.total_sample_size = int16(Result_test_covar.N);
S.model_select_index = int16(Result_test_covar.ind_for_model);
S.total_time = Result_test_covar.CPU_time_total;
disp(S)

end