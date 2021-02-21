close all;
clear all;

%% Parameter of the Gaussian function exp(-(a*x).^2)
a = 1;

%% 2D candidate grid on the square [lb_val, ub_val]^2
lb_val = 0;
ub_val = 1;

m = 50;
n = 50;
X = linspace(lb_val, ub_val, m);
Y = linspace(lb_val, ub_val, n);
XY = XY2D_grid(X,Y);
scatter(XY(:,1),XY(:,2));
xlim([lb_val, ub_val]);
ylim([lb_val, ub_val]);

pause(1);

%% Sequential Bayesian Quadrature (SBQ)
N = length(XY);
IIGau = SUB_GauK_db_int_2D_sq(a, lb_val, ub_val);
ZW = [];
WCE_arr = [];
Jmax = 50;
Jstp = 10;
for j=1:Jmax
    tmp_min = Inf;
    tmp_idx = 0;
    for i=1:(N-j+1)    
        tmp_ZW = [ZW; XY(i,:)];        
        [wc_en, w] = wce_ene(a, lb_val, ub_val, tmp_ZW);
        % [wc_en, w] = wce_ene_posw(a, lb_val, ub_val, tmp_ZW);
        if tmp_min > wc_en
            tmp_min = wc_en;
            tmp_idx = i;
        end
    end
    ZW = [ZW; XY(tmp_idx,:)];
    XY(tmp_idx,:) = [];
    
    if mod(j,Jstp)==0
        wce = IIGau + wc_en;
        WCE_arr = [WCE_arr; [j, wce]];
        display(wce);
    end
end

%% Output
display(w); % not necessarily positive!
display(wce);
scatter(ZW(:,1),ZW(:,2));
xlim([lb_val, ub_val]);
ylim([lb_val, ub_val]);

dlmwrite('Data_2D_SBQ_a1/data.txt', WCE_arr);
dlmwrite('Data_2D_SBQ_a1/data_nodes.txt', ZW);

%% Definitions of functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XY = XY2D_grid(X,Y)
    M = length(X);
    N = length(Y);
    XY = zeros(M*N,2);
    k = 0;
    for i=1:M
        for j=1:N
            k = k+1;
            XY(k,1) = X(i);
            XY(k,2) = Y(j);
        end
    end
end

function [wc_en, w] = wce_ene(a, lb_val, ub_val, XY)
    mat_K = exp(-a^2 * SUB_mat_dist2_2D(XY));
    vec_k = SUB_GauK_int_2D_sq(a, lb_val, ub_val, XY(:,1), XY(:,2));
    w = mat_K \ vec_k;
    wc_en = w' * mat_K * w - 2 * w' * vec_k;
end

function [wc_en, w] = wce_ene_posw(a, lb_val, ub_val, XY)
    mat_K = exp(-a^2 * SUB_mat_dist2_2D(XY));
    vec_k = SUB_GauK_int_2D_sq(a, lb_val, ub_val, XY(:,1), XY(:,2));

    N = length(XY(:,1));
    obj_func = @(w) (w' * mat_K * w - 2 * w' * vec_k);
    options = optimoptions(@fmincon, 'MaxFunctionEvaluations',1e+6);
    w = fmincon(obj_func, ones(N,1)/N, [],[],[],[], zeros(N,1)+0.005,[], [], options);

    wc_en = obj_func(w);        
end
