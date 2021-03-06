% 3D "mat_func" energy with a Gaussian external field for general weights w
function ene = SUB_ene_with_IntGauExtF_wght_3D(a, lb_val, ub_val, XY, w, mat_func)
    N = length(XY(:,1));
    mut_E = w' * mat_func(XY) * w;
    ext_F = - 2 * w' * SUB_GauK_int_3D_sq(a, lb_val, ub_val, XY(:,1), XY(:,2), XY(:,3));
    ene = mut_E + ext_F;
end
