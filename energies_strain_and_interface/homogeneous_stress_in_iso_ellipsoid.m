function sigma_i = homogeneous_stress_in_iso_ellipsoid( epsilon_c, epsilon_t, E, nu  )
% call: sigma_i = homogeneous_inclusion_stress( epsilon_c, epsilon_t, E, nu  )
% epsilon_c...coherent strain in inclusion - epsilon_c = S_ijkl epsilon_trans_ij

lame1 = nu*E / ((1+nu)*(1-2*nu));
G = E / 2*(1. + nu);

dev_ec = deviator( epsilon_c );
dev_et = deviator( epsilon_t );

% Equation 2.7 Eshelby 1957
sigma_i = lame1*(trace(epsilon_c) - trace(epsilon_t))*eye(3) + 2*G*(dev_ec - dev_et);

end

