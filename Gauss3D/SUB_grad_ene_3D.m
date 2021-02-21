% 2D gradient of the energy with respect to XY(i,:)
function grad_ene = SUB_grad_ene_3D(a, lb_val, ub_val, XY, idx, pdk_func)
    N = length(XY(:,1));
    x = XY(:,1);
    y = XY(:,2);
    z = XY(:,3);
    
    diff_x = zeros(N-1,1);
    diff_y = zeros(N-1,1);
    diff_z = zeros(N-1,1);
    dist2i = zeros(N-1,1);
    k = 1;
    for j=1:N
        if j~=idx
            diff_x(k) = x(idx)-x(j);
            diff_y(k) = y(idx)-y(j);            
            diff_z(k) = z(idx)-z(j);            
            dist2i(k) = diff_x(k)^2 + diff_y(k)^2 + diff_z(k)^2;
            k = k + 1;
        end
    end
    grad_Gau_x = 2 * sum(pdk_func(diff_x, dist2i));
    grad_Gau_y = 2 * sum(pdk_func(diff_y, dist2i));    
    grad_Gau_z = 2 * sum(pdk_func(diff_z, dist2i));    
    
    xx = x(idx);
    yy = y(idx);
    zz = z(idx);
    tmp_x = (erf(a*(ub_val-xx)) + erf(a*(xx-lb_val)));
    dtm_x = (2*a/sqrt(pi)) * ( -exp(-a^2 * (ub_val-xx)^2) + exp(-a^2 * (xx-lb_val)^2) );
    tmp_y = (erf(a*(ub_val-yy)) + erf(a*(yy-lb_val)));
    dtm_y = (2*a/sqrt(pi)) * ( -exp(-a^2 * (ub_val-yy)^2) + exp(-a^2 * (yy-lb_val)^2) );
    tmp_z = (erf(a*(ub_val-zz)) + erf(a*(zz-lb_val)));
    dtm_z = (2*a/sqrt(pi)) * ( -exp(-a^2 * (ub_val-zz)^2) + exp(-a^2 * (zz-lb_val)^2) );

    grad_ext_x = (sqrt(pi)/(2*a))^3 * dtm_x .* tmp_y .* tmp_z;
    grad_ext_y = (sqrt(pi)/(2*a))^3 * tmp_x .* dtm_y .* tmp_z;
    grad_ext_z = (sqrt(pi)/(2*a))^3 * tmp_x .* tmp_y .* dtm_z;
        
    grad_ene = (1/N^2) * [grad_Gau_x, grad_Gau_y, grad_Gau_z] - (2/N) * [grad_ext_x, grad_ext_y,  grad_ext_z];
end
