close all;
clear all;

%% Parameter of the Gaussian function exp(-(a*x).^2)
a = 1;

%% 3D initial grid on the square [lb_val, ub_val]^2
lb_val = 0;
ub_val = 1;
ini_mrgn = 0.1; 

N_arr = [10:10:100];
lngth = length(N_arr);
wce_arr = zeros(1,lngth);

for ell = 1:length(N_arr)
    N = N_arr(ell);
    display(N);
    XYZ = (lb_val + ini_mrgn) + (ub_val - lb_val - 2 * ini_mrgn) * rand(N,3);

    scatter3(XYZ(:,1),XYZ(:,2), XYZ(:,3));
    w = ones(N,1)/N;
    display(SUB_ene_with_IntGauExtF_wght_3D(a, lb_val, ub_val, XYZ, w, @(xy)mat_Gau_3D(a,xy)));
    display(SUB_ene_with_IntGauExtF_wght_3D(a, lb_val, ub_val, XYZ, w, @mat_m2rd_3D));

    pause(1);

    %% Point-wise gradient descent of an energy
    default_gamma = 1;
    
    pwr = 1.25; 
    mrgn = 0.12;
    % mrgn = log10(N)/10;
    
    norm_grad_eps = 1e-4;
    max_itr = 10000;
    for k = 1:max_itr
        max_norm_grad = 0;
        for idx=1:N

        %     grad_ene = SUB_grad_ene_3D(a, lb_val, ub_val, XYZ, idx, @(s,t)pdk_Gau(a,s,t)); % Gaussian
            grad_ene = SUB_grad_ene_3D(a, lb_val, ub_val, XYZ, idx, @pdk_rd); % Logarithmic

        %     grad_reg = [0 0 0];
            grad_reg = (1/N^pwr) * SUB_reg_grad_rec_dist_3D_sq(XYZ(idx,1), XYZ(idx,2), XYZ(idx,3), lb_val-mrgn, ub_val+mrgn);

            grad = grad_ene + grad_reg;

            gamma = SUB_set_step_size_3D_sq(XYZ(idx,:), grad, default_gamma, lb_val, ub_val);    
            XYZ(idx,:) = XYZ(idx,:) - gamma * grad;

            if max_norm_grad < norm(grad)
                max_norm_grad = norm(grad);
            end        
        end
        if max_norm_grad < norm_grad_eps
            break;
        end
    end

    scatter3(XYZ(:,1),XYZ(:,2), XYZ(:,3));
    display(SUB_ene_with_IntGauExtF_wght_3D(a, lb_val, ub_val, XYZ, w, @(xyz)mat_Gau_3D(a,xyz)));
    display(SUB_ene_with_IntGauExtF_wght_3D(a, lb_val, ub_val, XYZ, w, @mat_m2rd_3D));

    pause(1);

    %% Optimal weights
    mat_K = exp(-a^2 * SUB_mat_dist2_3D(XYZ));
    vec_k = SUB_GauK_int_3D_sq(a, lb_val, ub_val, XYZ(:,1), XYZ(:,2), XYZ(:,3));

    w = mat_K \ vec_k;

    fin_Gau_ene = SUB_ene_with_IntGauExtF_wght_3D(a, lb_val, ub_val, XYZ, w, @(xyz)mat_Gau_3D(a,xyz));
    display(fin_Gau_ene);

    IIGau = SUB_GauK_db_int_3D_sq(a, lb_val, ub_val);
    fin_WCE = IIGau + fin_Gau_ene;
    display(fin_WCE);
    
    %% Output
    wce_arr(ell) = fin_WCE;
end

dlmwrite('Data_3D_PWGD_a1/data.txt', [N_arr', wce_arr']);
dlmwrite('Data_3D_PWGD_a1/data_nodes.txt', XYZ);

%% Definitions of functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [[ Gaussian ]]

% 3D kernel matrix of the Gaussian function
function mG_3D = mat_Gau_3D(a, XYZ)
    dist2 = SUB_mat_dist2_3D(XYZ);
    mG_3D = exp(-a^2 * dist2);
end

% Partial derivative of the Gaussian kernel
function pdxy = pdk_Gau(a, diff_xyz, dist2i)
    pdxy = -2 * a^2 .* diff_xyz .* exp(-a^2 * dist2i); 
end

% [[ Logarithmic ]]

% 3D kernel matrix of the reciprocals of the distances
function m2rd_3D = mat_m2rd_3D(XYZ)
    N = length(XYZ(:,1));
    dist2 = SUB_mat_dist2_3D(XYZ);
    for i=1:N
        dist2(i,i) = 1;
    end    
    
    m2rd_3D = 2 * sqrt(pi) ./ sqrt(dist2);
    for i=1:N
        m2rd_3D(i,i) = 0;
    end    
end

% Partial derivative of the reciprocal distance kernel
function pdxyz = pdk_rd(diff_xyz, dist2i)
    pdxyz = - 2 * sqrt(pi) * diff_xyz ./ ((dist2i).^(3/2));
end
