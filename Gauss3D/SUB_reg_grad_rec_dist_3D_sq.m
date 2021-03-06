% Differential of a regularization term by a reciprocal of the distance functions
function grad_xy = SUB_reg_grad_rec_dist_3D_sq(x, y, z, lb_val, ub_val)
    
    tmp_x = diff_psi(x, lb_val, ub_val);
    tmp_y = diff_psi(y, lb_val, ub_val);    
    tmp_z = diff_psi(z, lb_val, ub_val);    
        
    grad_xy = [tmp_x, tmp_y, tmp_z];
end

function z = diff_psi(x, lb_val, ub_val)
    z = -1./(x-lb_val).^2 + 1./(ub_val-x).^2;
end

function z = psi(x, lb_val, ub_val)
    z = 1./(x-lb_val) + 1./(ub_val-x);
end
