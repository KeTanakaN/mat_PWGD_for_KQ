% Rule for a step size such that p is in [lb_val, ub_val]^3
function gamma = SUB_set_step_size_3D_sq(p, grad, default_gamma, lb_val, ub_val)
    ub_gamma = default_gamma * ones(1,3);
    for i=1:3
        if grad(i) > 0
            ub_gamma(i) = (p(i) - lb_val) / grad(i);
        elseif grad(i) < 0
            ub_gamma(i) = (p(i) - ub_val) / grad(i);
        end
    end    
    
    gamma = min([default_gamma, ub_gamma]);
end
