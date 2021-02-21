% Differential of a regularization term by a product of logarithmic functions
function grad_xy = SUB_reg_grad_log_2D_sq(x, y, lb_val, ub_val)

%     tmp_x = diff_psi(x, lb_val, ub_val) .* psi(x, lb_val, ub_val);
%     tmp_y = psi(y, lb_val, ub_val) .*  diff_psi(y, lb_val, ub_val);    
    
    tmp_x = diff_psi(x, lb_val, ub_val);
    tmp_y = diff_psi(y, lb_val, ub_val);    
        
    grad_xy = [tmp_x, tmp_y];
end

function z = diff_psi(x, lb_val, ub_val)
%     tmp = -1./(x-lb_val) + 1./(ub_val-x);
%     z = tmp .* psi(x, lb_val, ub_val);

    z = -1./(x-lb_val) + 1./(ub_val-x);
end

function z = psi(x, lb_val, ub_val)
    z = -log((x-lb_val).*(ub_val-x));
end
