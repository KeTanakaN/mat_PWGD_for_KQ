close all;
clear all;

%% Parameter of the Gaussian function exp(-(a*x).^2)
a = 1;

%% 2D initial grid on the square [lb_val, ub_val]^2
lb_val = 0;
ub_val = 1;
ini_mrgn = 0.1; 

N_arr = [10:10:50];
lngth = length(N_arr);
wce_arr = zeros(1,lngth);

for ell = 1:length(N_arr)
    N = N_arr(ell);
    display(N);
    XY = (lb_val + ini_mrgn) + (ub_val - lb_val - 2 * ini_mrgn) * rand(N,2);

    scatter(XY(:,1),XY(:,2));
    w = ones(N,1)/N;
    display(SUB_ene_with_IntGauExtF_wght_2D(a, lb_val, ub_val, XY, w, @(xy)mat_Gau_2D(a,xy)));
    display(SUB_ene_with_IntGauExtF_wght_2D(a, lb_val, ub_val, XY, w, @mat_m2log_2D));

    % pause(1);

    %% Point-wise gradient descent of an energy
    default_gamma = 0.1;
    
    pwr = 0.50; 
    mrgn = 0.50;
    % mrgn = log10(N)/10;
    
    norm_grad_eps = 1e-4;
    max_itr = 1000;
    for k = 1:max_itr
        max_norm_grad = 0;
        for idx=1:N

    %         grad_ene = SUB_grad_ene_2D(a, lb_val, ub_val, XY, idx, @(s,t)pdk_Gau(a,s,t)); % Gaussian
            grad_ene = SUB_grad_ene_2D(a, lb_val, ub_val, XY, idx, @pdk_log); % Logarithmic

    %         grad_reg = [0 0];
            grad_reg = (1/N^pwr) * SUB_reg_grad_log_2D_sq(XY(idx,1), XY(idx,2), lb_val-mrgn, ub_val+mrgn);

            grad = grad_ene + grad_reg;

            gamma = SUB_set_step_size_2D_sq(XY(idx,:), grad, default_gamma, lb_val, ub_val);    
            XY(idx,:) = XY(idx,:) - gamma * grad;

            if max_norm_grad < norm(grad)
                max_norm_grad = norm(grad);
            end        
        end
        if max_norm_grad < norm_grad_eps
            break;
        end
    end

    scatter(XY(:,1),XY(:,2));
    display(SUB_ene_with_IntGauExtF_wght_2D(a, lb_val, ub_val, XY, w, @(xy)mat_Gau_2D(a,xy)));
    display(SUB_ene_with_IntGauExtF_wght_2D(a, lb_val, ub_val, XY, w, @mat_m2log_2D));

    pause(1);

    %% Optimal weights
    mat_K = exp(-a^2 * SUB_mat_dist2_2D(XY));
    vec_k = SUB_GauK_int_2D_sq(a, lb_val, ub_val, XY(:,1), XY(:,2));

    w = mat_K \ vec_k; % negative weights appear -> constrained optimization is needed

%     obj_func = @(w) (w' * mat_K * w - 2 * w' * vec_k);
%     options = optimoptions(@fmincon, 'MaxFunctionEvaluations',1e+6);
%     w = fmincon(obj_func, ones(N,1)/N, [],[],[],[], zeros(N,1),[], [], options);

    fin_Gau_ene = SUB_ene_with_IntGauExtF_wght_2D(a, lb_val, ub_val, XY, w, @(xy)mat_Gau_2D(a,xy));
    display(fin_Gau_ene);

    IIGau = SUB_GauK_db_int_2D_sq(a, lb_val, ub_val);
    fin_WCE = IIGau + fin_Gau_ene;
    display(fin_WCE);
    
    %% Output
    wce_arr(ell) = fin_WCE;
end

dlmwrite('Data_2D_PWGD_a1/data.txt', [N_arr', wce_arr']);
dlmwrite('Data_2D_PWGD_a1/data_nodes.txt', XY);

%% Definitions of functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [[ Gaussian ]]

% 2D kernel matrix of the Gaussian function
function mG_2D = mat_Gau_2D(a, XY)
    dist2 = SUB_mat_dist2_2D(XY);
    mG_2D = exp(-a^2 * dist2);
end

% Partial derivative of the Gaussian kernel
function pdxy = pdk_Gau(a, diff_xy, dist2i)
    pdxy = -2 * a^2 .* diff_xy .* exp(-a^2 * dist2i); 
end

% [[ Logarithmic ]]

% 2D kernel matrix of the logarithm
function m2l_2D = mat_m2log_2D(XY)
    N = length(XY(:,1));
    dist2 = SUB_mat_dist2_2D(XY);
    for i=1:N
        dist2(i,i) = 1;
    end    
    m2l_2D = -4 * (log(dist2)/2);
end

% Partial derivative of the logarithmic kernel
function pdxy = pdk_log(diff_xy, dist2i)
    pdxy = -4 * diff_xy ./ dist2i;
end
