clear all
close all

addpath('./Other_functions')
addpath('./spinterp_v5.1.1_vec','./spinterp_v5.1.1_vec/private_function')
addpath('./Toy_Surrogates')
%% ============= Parameters ===============================================
% Parameters
perturb_ratio1 = 0.03;
perturb_ratio2 = 0.03;
d=2;
range = [[-1,1]*perturb_ratio1;[-1,1]*perturb_ratio2];

tol = 1e-2; % tolerance for stopping
interpolate=false;
norm_l = 'L2err_toy';
Eps=1e-3;
theta=0.5;





%--------------------------------------------------------------------------
%load hfm mesh and common mesh
hfm_mesh_level = 8;
com_mesh_level = 8;
index_k = 2:hfm_mesh_level;
NoLFM = length(index_k);

[Mesh_com,Mesh_high,h_high,QoI_h] = load_hfm_mesh_n_com_mesh(com_mesh_level,hfm_mesh_level);


mean_hfm = zeros(QoI_h,NoLFM);
var_hfm = zeros(1,NoLFM);
mean_hfm_temp =0;
var_hfm_temp = 0;
covar = zeros(1,NoLFM);

mean_lfm = zeros(QoI_h,NoLFM);
var_lfm = zeros(1,NoLFM);
var_lfm_temp = zeros(1,NoLFM);
QoI_len = zeros(1,NoLFM);



% pre-sample for mean: both hf and lf surrog.
NN = 3;
rand_sample_1 = [(rand(NN,1)*2-1)*perturb_ratio1,(rand(NN,1)*2-1)*perturb_ratio2];

for i=1:NN
    u = FEM_solver(hfm_mesh_level,rand_sample_1(i,1),rand_sample_1(i,2));
    mean_hfm_temp = mean_hfm_temp+u;
    var_hfm_temp = var_hfm_temp + feval(norm_l,u,0,h_high)^2;
end
mean_hfm_temp = mean_hfm_temp/NN;
mean_hfm =repmat(mean_hfm_temp,1,NoLFM);
norm_scale = feval(norm_l,u,0,h_high);

if interpolate
    % mean_hfm = interp2grid(mean_hfm,Mesh_high,Mesh_com);
end



for i=1:NN
    for k=index_k
        %load lfm mesh and surrogate
        [Mesh_low,QoI_low,h_low,z] = load_lfm(k,interpolate);
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



%----------------------------------------------------------------------
% Dynamic sampling to find statitics
%----------------------------------------------------------------------
sum_time=zeros(1,length(index_k)+1);
sum_time_proj = zeros(1,8);
sum_time_proj_fine = 0;
t_total = tic;

dN=10; % batch of sample
Addsample=true;
C = zeros(1,length(index_k)+1);
sample_stat_vec = zeros(10, 4+length(C),length(index_k));
p=0;


while Addsample

    for k=index_k

        for i =1:dN

            j=p+i;

            % generate sample
            rand_sample = [(rand(1,1)*2-1)*perturb_ratio1,(rand(1,1)*2-1)*perturb_ratio2];

            %----------------------------------------------------------------------
            %hfm--soln of discretized pde
            t1 = tic;
            u = FEM_solver(hfm_mesh_level,rand_sample(1),rand_sample(2));
            tt = toc(t1);
            sum_time(1) = sum_time(1)+tt;
            %mean
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
            var_hfm(k-1) = var_hfm(k-1) + feval(norm_l,d_mean_hfm, u-mean_hfm(:,k-1),h_high); % var
            %----------------------------------------------------------------------


            %----------------------------------------------------------------------
            %lfm--soln of surrogate
            t3 = tic;
            u_surrog = Surrog_Eval(z_store{k-1},range,d,QoI_len(k-1),rand_sample);
            tt = toc(t3);
            sum_time(k) = sum_time(k)+tt;

            %mean
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



            % store 2 latest parameters into sample_stat_vec
            if i>8

                C = sum_time/j;
                sigma1 = sqrt(var_hfm(k-1)/(j-1));
                sigmak = sqrt(var_lfm(k-1)/(j-1));
                rho = (covar(k-1)/(j-1))./(sqrt(var_hfm(k-1)/(j-1))*sqrt(var_lfm(k-1)/(j-1)));
                sample_stat_vec(i-8,:,k-1) = [sigma1, sigmak, rho, covar(k-1)/(j-1),C];  %

            end

        end
    end




    %----------------------------------------------------------------------
    % model selection
    %----------------------------------------------------------------------
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
    V_MF_j_1 = mean(sigma1)^2/N_j(1)+ sum((1./N_j(1:end-1) - 1./N_j(2:end)).*(alpha_j.^2.*sigmak(2:end).^2-2*alpha_j.*rho(2:end).*mean(sigma1).*sigmak(2:end)));




    % stopping criterion for sufficient sample size
    relative_error = abs(xi_j - xi_j_1)/abs(xi_j); % if sigma_1, rho_k satisfy the convergence

    if relative_error < tol
        Addsample = false;
    else
        Addsample = true;
        p=p+dN;
    end


end
t_total = toc(t_total);

N=j;
Result_test_covar.N = N;
Result_test_covar.var_hfm = sqrt(var_hfm/(N-1)); % sigma1--standdard deviation of hfm
Result_test_covar.var_lfm = sqrt(var_lfm/(N-1)); % sigmak--standdard deviation of lfm
Result_test_covar.covar = covar/(N-1); % covariance btw hfm & lfm
Result_test_covar.rho = (covar/(N-1))./(sqrt(var_hfm/(N-1)).*sqrt(var_lfm/(N-1))); % rho
Result_test_covar.CPU_time_fine_direct = sum_time/N;
Result_test_covar.CPU_time_total = t_total;
Result_test_covar.sum_time_proj = sum_time_proj;
Result_test_covar.sum_time_proj_fine = sum_time_proj_fine;
Result_test_covar.sample_stat_vec = sample_stat_vec;

