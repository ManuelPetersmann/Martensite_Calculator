function [res_B, res_vec, RS] = bain_opt_shear( B, dS, S_accummulated, u )
% try to find a shear that is similar to the bain strain - I dont think
% this is a good idea anymore ^^

[US, RS] = polardecomposition( dS * S_accummulated );

res_B = sum( sum( abs(B - US )));

u2_dS = dS * S_accummulated * u; 
new_R = rotation_between_vectors( u2_dS, u );
u2_dS_unrot = new_R * u2_dS;

res_vec = norm( u - u2_dS_unrot ); % here the euklidean norm is taken! - c.f. leastsquares problem

LT    = new_R * martensite.U;
[ theta_cp_new1, ~ ] = min_misorientation( austenite.CPPs, LT1, true );
end

