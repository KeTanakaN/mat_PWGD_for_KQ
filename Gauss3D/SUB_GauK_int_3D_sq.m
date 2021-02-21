% 3D Integral of the Gaussian kernel on [lb_val, ub_val]^2
function z = SUB_GauK_int_3D_sq(a, lb_val, ub_val, x, y, z)
    tmp_x = (erf(a*(ub_val-x)) + erf(a*(x-lb_val)));
    tmp_y = (erf(a*(ub_val-y)) + erf(a*(y-lb_val)));
    tmp_z = (erf(a*(ub_val-z)) + erf(a*(z-lb_val)));
    z = (sqrt(pi)/(2*a))^3 * tmp_x .* tmp_y .* tmp_z; 
end
