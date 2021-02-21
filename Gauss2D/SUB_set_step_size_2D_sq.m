% Rule for a step size such that p is in [lb_val, ub_val]^2
function gamma = SUB_set_step_size_2D_sq(p, grad, default_gamma, lb_val, ub_val)
    ub_gamma = default_gamma * ones(1,2);
    for i=1:2
        if grad(i) > 0
            ub_gamma(i) = (p(i) - lb_val) / grad(i);
        elseif grad(i) < 0
            ub_gamma(i) = (p(i) - ub_val) / grad(i);
        end
    end    
    
    gamma = min([default_gamma, ub_gamma]);
end
