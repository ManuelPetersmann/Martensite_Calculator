function [new_res, new_R] = perp_ILS( US, dS, u)
u2_dS = US * dS * u;
new_R = rotation_between_vectors( u2_dS, u );
u2_dS_unrot = new_R * u2_dS;
new_res = norm( u - u2_dS_unrot );
end
