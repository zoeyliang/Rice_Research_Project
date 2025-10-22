clear all
close all
format shortE

% p = genpath('/Users/jiaxingliang/Dropbox/FEEQS.M-0.1p/');
% p = genpath('/Toy');
% addpath(p)


addpath('./Other_functions')
addpath('./spinterp_v5.1.1_vec','./spinterp_v5.1.1_vec/private_function')
addpath('./Toy_Surrogates')


%% ============= Parameters ===============================================
% Parameters
perturb_ratio1 = 0.03; % percentage of the first parameter.
perturb_ratio2 = 0.03; % percentage of the second parameter.
d=2; % dimension of the parameter space
range = [[-1,1]*perturb_ratio1;[-1,1]*perturb_ratio2]; % range of parameter space, size dx1.

tol = 1e-2; % tolerance for stopping
interpolate = false; %whether the surrogate has already been interpolated to the common mesh or not.
norm_l = 'L2_inprod_toy'; %inner product function
Eps=1e-3; % threshold for nMSE.
theta=0.5; % splitting ratio between discretization and statistical error.





%--------------------------------------------------------------------------
%load hfm mesh and common mesh
hfm_mesh_level = 8; % level of high fidelity mesh
com_mesh_level = 8; % level of common mesh
index_k = 2:8; % index for low fidelity models
NoLFM = length(index_k); % # of low fidelity models
[Mesh_com,Mesh_high,h_high,QoI_h] = load_hfm_mesh_n_com_mesh_toy(com_mesh_level,hfm_mesh_level);


%initialization
mean_hfm = zeros(QoI_h,NoLFM);
var_hfm = zeros(1,NoLFM);
mean_hfm_temp =0;
var_hfm_temp = 0;
covar = zeros(1,NoLFM);

mean_lfm = zeros(QoI_h,NoLFM);
var_lfm = zeros(1,NoLFM);
var_lfm_temp = zeros(1,NoLFM);
QoI_len = zeros(1,NoLFM);


%--------------------------------------------------------------------------
% pre-sample for mean: both hf and lf surrog.
%--------------------------------------------------------------------------
NN = 3;
rand_sample_1 = [(rand(NN,1)*2-1)*perturb_ratio1,(rand(NN,1)*2-1)*perturb_ratio2];

for i=1:NN
    u = FEM_solver(hfm_mesh_level,rand_sample_1(i,1),rand_sample_1(i,2));
    mean_hfm_temp = mean_hfm_temp+u;
    var_hfm_temp = var_hfm_temp + feval(norm_l,u,u,h_high)^2;
end
mean_hfm_temp = mean_hfm_temp/NN;
mean_hfm =repmat(mean_hfm_temp,1,NoLFM);
norm_scale = sqrt(feval(norm_l,mean_hfm_temp,mean_hfm_temp,h_high));

% this is only useful when the low fidelity model on coarser grid is
% interpolated to the common fine grid as the high fidelity model.
if interpolate
    % mean_hfm = interp2grid(mean_hfm,Mesh_high,Mesh_com);
end

for i=1:NN
    for k=index_k
        %load lfm mesh and surrogate
        [Mesh_low,QoI_low,h_low,z] = load_lfm_toy(k,interpolate);
        z=z.z;
        z_store{k-1} = z;
        Mesh_low_store{k-1} = Mesh_low;
        h_low_store{k-1} = h_low;
        % QoI_low_store{k} = QoI_low;


        u_surrog = Surrog_Eval(z,range,d,QoI_low,rand_sample_1(i,:));
        QoI_len(k-1) = size(u_surrog,1);
        mean_lfm(1:QoI_len(k-1),k-1) = mean_lfm(1:QoI_len(k-1),k-1)+u_surrog;
        var_lfm_temp(:,k-1) = var_lfm_temp(:,k-1) + feval(norm_l,u_surrog,u_surrog,h_low)^2;
    end
end
mean_lfm = mean_lfm/NN;



%--------------------------------------------------------------------------
% Dynamic sampling to find statitics
%--------------------------------------------------------------------------
sum_time=zeros(1,NoLFM+1);
sum_time_proj = zeros(1,NoLFM+1);
sum_time_proj_fine = 0;
C = zeros(1,NoLFM+1);
sample_stat_vec = zeros(10, 4+length(C),NoLFM);

dN=10; % batch of sample
Addsample=true; % flag to decide whether to add samples or not.
p=0; % index for batch of samples



t_total = tic;
while Addsample

    for k=index_k % loop over all low fidelity models

        for i =1:dN % loop over a batch of samples

            j=p+i;

            % generate sample
            rand_sample = [(rand(1,1)*2-1)*perturb_ratio1,(rand(1,1)*2-1)*perturb_ratio2];

            %----------------------------------------------------------------------
            %hfm--soln of discretized pde
            t1 = tic;
            u = FEM_solver(hfm_mesh_level,rand_sample(1),rand_sample(2));
            tt = toc(t1);
            sum_time(1) = sum_time(1)+tt;


            if interpolate
                % t2 = tic;
                % u = interp2grid(u,Mesh_high,Mesh_com);
                % tt2 = toc(t2);
                % sum_time_proj_fine = sum_time_proj_fine+tt2;
            end

            % Welford update statistics
            d_mean_hfm = u - mean_hfm(:,k-1);
            mean_hfm(:,k-1) = mean_hfm(:,k-1) + d_mean_hfm/j; % mean

            % values to generate var and correlation
            var_hfm(k-1) = var_hfm(k-1) + feval(norm_l,d_mean_hfm, u-mean_hfm(:,k-1),h_high); % variance
            %----------------------------------------------------------------------


            %----------------------------------------------------------------------
            %lfm--soln of surrogate
            t3 = tic;
            u_surrog = Surrog_Eval(z_store{k-1},range,d,QoI_len(k-1),rand_sample);
            tt = toc(t3);
            sum_time(k) = sum_time(k)+tt;

            % Welford update statistics
            d_mean_lfm = u_surrog(1:QoI_len(k-1)) - mean_lfm(1:QoI_len(k-1),k-1);
            mean_lfm(1:QoI_len(k-1),k-1) = mean_lfm(1:QoI_len(k-1),k-1) + d_mean_lfm/j; % mean
            % values to generate var and correlation
            var_lfm(k-1) = var_lfm(k-1) + feval(norm_l,d_mean_lfm, u_surrog-mean_lfm(1:QoI_len(k-1),k-1),h_low_store{k-1}); % var
            %----------------------------------------------------------------------


            % covariance btw hfm & lfm
            if interpolate
                % d_mean_hfm_interp = d_mean_hfm(1:QoI);
                % d_mean_lfm_interp = u_surrog(1:QoI)-mean_lfm(1:QoI,k-1)
            else

                t2 = tic;
                d_mean_hfm_interp = d_mean_hfm; %interp2grid(d_mean_hfm(1:QoI),Mesh_high,Mesh_com);
                d_mean_lfm_interp = interp2grid_toy(u_surrog-mean_lfm(1:QoI_len(k-1),k-1),Mesh_low_store{k-1},Mesh_com);
                tt2 = toc(t2);
                sum_time_proj(k-1) = sum_time_proj(k-1)+tt2;
            end

            covar(k-1) = covar(k-1) + feval(norm_l,d_mean_hfm_interp,d_mean_lfm_interp,h_high);



            % store 2 or more latest parameters into sample_stat_vec
            if i>dN-2

                C = sum_time/j;
                sigma1 = sqrt(var_hfm(k-1)/(j-1));
                sigmak = sqrt(var_lfm(k-1)/(j-1));
                rho = (covar(k-1)/(j-1))./(sqrt(var_hfm(k-1)/(j-1))*sqrt(var_lfm(k-1)/(j-1)));
                sample_stat_vec(i-8,:,k-1) = [sigma1, sigmak, rho, covar(k-1)/(j-1),C];
            end

        end
    end




    %----------------------------------------------------------------------
    % model selection
    %----------------------------------------------------------------------
    % second latest update
    sigma1 = squeeze(sample_stat_vec(1,1,:))';
    sigmak = [mean(sigma1);squeeze(sample_stat_vec(1,2,:))]';
    rho = [1;squeeze(sample_stat_vec(1,3,:))]';

    C=squeeze(sample_stat_vec(1,5:end,end));

    [ind_for_model, sum_sqrt_C_rho] = MFMC_model_selection_exhausted(rho,C);
    rho = rho(ind_for_model);
    rho_p1 = [rho(2:end),0];
    sigmak = sigmak(ind_for_model);
    C = C(ind_for_model);
    % compute xi^{j-1}/alpha_{j-1}/N_{j-1}/V_MF_{j-1} by (21) with the selected K^* models
    xi_j_1 = 1/C(1)*sum(sqrt((rho.^2- rho_p1.^2).*C))^2;
    alpha_j_1 = rho(2:end).*mean(sigma1)./sigmak(2:end);
    N_j_1 = max(ceil(mean(sigma1)^2/(norm_scale^2*Eps^2*(1-theta))*sqrt((rho.^2- rho_p1.^2)./C)*sum_sqrt_C_rho),2);
    V_MF_j_1 = mean(sigma1)^2/N_j_1(1)+ sum((1./N_j_1(1:end-1) - 1./N_j_1(2:end)).*(alpha_j_1.^2.*sigmak(2:end).^2-2*alpha_j_1.*rho(2:end).*mean(sigma1).*sigmak(2:end)));
    




    %latest update
    sigma1 = squeeze(sample_stat_vec(2,1,:))';
    sigmak = [mean(sigma1);squeeze(sample_stat_vec(2,2,:))]';
    rho = [1;squeeze(sample_stat_vec(2,3,:))]';

    C=squeeze(sample_stat_vec(2,5:end,end));

    [ind_for_model, sum_sqrt_C_rho] = MFMC_model_selection_exhausted(rho,C);
    rho = rho(ind_for_model);
    rho_p1 = [rho(2:end),0];
    sigmak = sigmak(ind_for_model);
    C = C(ind_for_model);
    % compute efficienty xi^{j}/ weights alpha_j/ sample size N_j/ MFMC variance V_MF_j by (21) with the selected K^* models
    xi_j = 1/C(1)*sum(sqrt((rho.^2- rho_p1.^2).*C))^2;
    alpha_j = rho(2:end).*mean(sigma1)./sigmak(2:end);
    N_j = max(ceil(mean(sigma1)^2/(norm_scale^2*Eps^2*(1-theta))*sqrt((rho.^2- rho_p1.^2)./C)*sum_sqrt_C_rho),2);
    V_MF_j = mean(sigma1)^2/N_j(1)+ sum((1./N_j(1:end-1) - 1./N_j(2:end)).*(alpha_j.^2.*sigmak(2:end).^2-2*alpha_j.*rho(2:end).*mean(sigma1).*sigmak(2:end)));
    %Note: the true V_MF=theta*Eps^2*norm_scale^2.
    %----------------------------------------------------------------------

    % stopping criterion for sufficient sample size
    relative_error = abs(xi_j - xi_j_1)/abs(xi_j); % if sigma_1, rho_k satisfy the convergence


    % decide if extra samples are needed or not
    if relative_error < tol
        Addsample = false;
    else
        Addsample = true;
        p=p+dN; % add a batch of dN samples
    end



end
t_total = toc(t_total);





% Before selection---------------------------------------------------------
Result_test_covar.BeforeSelection_sigmak = [mean(sigma1);squeeze(sample_stat_vec(2,2,:))]';
Result_test_covar.BeforeSelection_rho = [1;squeeze(sample_stat_vec(2,3,:))]';
Result_test_covar.BeforeSelection_C = squeeze(sample_stat_vec(2,5:end,end));
Result_test_covar.BeforeSelection_covar = covar/(j-1); % covariance btw hfm & lfm

% After selection----------------------------------------------------------
Result_test_covar.ind_for_model = ind_for_model;
Result_test_covar.sigmak = sigmak; % sigma1--standdard deviation of hfm % sigmak--standdard deviation of lfm
Result_test_covar.rho = rho; % rho
Result_test_covar.C = C;
%latest update
Result_test_covar.xi_j = xi_j; % MF to MC efficiency
Result_test_covar.alpha_j = alpha_j; % coefficients in the MFMC estimator for correction.
Result_test_covar.N_j = N_j; % sample size
Result_test_covar.V_MF_j = V_MF_j; % estimated V_MF
Result_test_covar.V_MF = theta*Eps^2*norm_scale^2; % true V_MF
% second latest update
Result_test_covar.xi_j_1 = xi_j_1 ;
Result_test_covar.alpha_j_1 = alpha_j_1;
Result_test_covar.N_j_1 = N_j_1;
Result_test_covar.V_MF_j_1 = V_MF_j_1;


Result_test_covar.N = j; % total sample size
Result_test_covar.CPU_time_fine_direct = sum_time/j;
Result_test_covar.CPU_time_total = t_total;
Result_test_covar.sum_time_proj = sum_time_proj;
Result_test_covar.sum_time_proj_fine = sum_time_proj_fine;
Result_test_covar.sample_stat_vec = sample_stat_vec;

plot_n_print(Result_test_covar)







