function [new_res, new_R] = perp_ILS( SU, dS,  u)
%% assumption - on the incoherent interface the Bain strain evolves and the stresses are immediately relaxed by the surface dislocations!
u2_dS = dS * SU * u; 
% I also had: shear must be added from right since we solve it as inverse problem.
%% actually this is stupid - the rotation could get bigger and it would look like the shear achieves the improvement
new_R = rotation_between_vectors( u2_dS, u );
%u2_dS_unrot = new_R * u2_dS;
%new_res = norm( u - u2_dS_unrot ); % here the euklidean norm is taken! - c.f. leastsquares problem
new_res = norm( u - u2_dS );
end
