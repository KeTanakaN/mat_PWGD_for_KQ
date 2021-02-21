close all;
clear all;

%% Parameter of the Gaussian function exp(-(a*x).^2)
a = 1;

%% 3D candidate grid on the square [lb_val, ub_val]^2
lb_val = 0;
ub_val = 1;

l = 20;
m = 20;
n = 20;
X = linspace(lb_val, ub_val, l);
Y = linspace(lb_val, ub_val, n);
Z = linspace(lb_val, ub_val, m);
XYZ = XY3D_grid(X,Y,Z);
scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3));
xlim([lb_val, ub_val]);
ylim([lb_val, ub_val]);
zlim([lb_val, ub_val]);

pause(1);

%% Sequential Bayesian Quadrature (SBQ)
N = length(XYZ);
IIGau = SUB_GauK_db_int_3D_sq(a, lb_val, ub_val);
ZW = [];
WCE_arr = [];
Jmax = 100;
Jstp = 10;
for j=1:Jmax
    tmp_min = Inf;
    tmp_idx = 0;
    for i=1:(N-j+1)    
        tmp_ZW = [ZW; XYZ(i,:)];        
        [wc_en, w] = wce_ene(a, lb_val, ub_val, tmp_ZW);
        % [wc_en, w] = wce_ene_posw(a, lb_val, ub_val, tmp_ZW);
        if tmp_min > wc_en
            tmp_min = wc_en;
            tmp_idx = i;
        end
    end
    ZW = [ZW; XYZ(tmp_idx,:)];
    XYZ(tmp_idx,:) = [];
    
    if mod(j,Jstp)==0
        wce = IIGau + wc_en;
        WCE_arr = [WCE_arr; [j, wce]];
        display(wce);
    end
end

%% Output
display(w); % not necessarily positive!
display(wce);
scatter3(ZW(:,1),ZW(:,2),ZW(:,3));
xlim([lb_val, ub_val]);
ylim([lb_val, ub_val]);
zlim([lb_val, ub_val]);

dlmwrite('Data_3D_SBQ_a1/data.txt', WCE_arr);
dlmwrite('Data_3D_SBQ_a1/data_nodes.txt', ZW);

%% Definitions of functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XYZ = XY3D_grid(X,Y,Z)
    L = length(X);
    M = length(Y);
    N = length(Z);
    XYZ = zeros(L*M*N,3);
    k = 0;
    for l=1:L
        for m=1:M
            for n=1:N
                k = k+1;
                XYZ(k,1) = X(l);
                XYZ(k,2) = Y(m);
                XYZ(k,3) = Z(n);
            end
        end
    end
end

function [wc_en, w] = wce_ene(a, lb_val, ub_val, XYZ)
    mat_K = exp(-a^2 * SUB_mat_dist2_3D(XYZ));
    vec_k = SUB_GauK_int_3D_sq(a, lb_val, ub_val, XYZ(:,1), XYZ(:,2), XYZ(:,3));
    w = mat_K \ vec_k;
    wc_en = w' * mat_K * w - 2 * w' * vec_k;
end

function [wc_en, w] = wce_ene_posw(a, lb_val, ub_val, XYZ)
    mat_K = exp(-a^2 * SUB_mat_dist2_3D(XYZ));
    vec_k = SUB_GauK_int_3D_sq(a, lb_val, ub_val, XYZ(:,1), XYZ(:,2), XYZ(:,3));

    N = length(XYZ(:,1));
    obj_func = @(w) (w' * mat_K * w - 2 * w' * vec_k);
    options = optimoptions(@fmincon, 'MaxFunctionEvaluations',1e+6);
    w = fmincon(obj_func, ones(N,1)/N, [],[],[],[], zeros(N,1)+0.005,[], [], options);

    wc_en = obj_func(w);        
end
