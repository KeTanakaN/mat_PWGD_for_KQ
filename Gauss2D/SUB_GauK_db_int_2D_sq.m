% 2D double integral of the Gaussian kernel on [lb_val, ub_val]^2
function z = SUB_GauK_db_int_2D_sq(a, lb_val, ub_val)
    tmp_1 = exp(-(a * (lb_val - ub_val))^2);
    tmp_2 = sqrt(pi) * a * (lb_val - ub_val) * erf(a * (lb_val - ub_val));
    z = ((tmp_1 + tmp_2 - 1)/a^2)^2;
end
