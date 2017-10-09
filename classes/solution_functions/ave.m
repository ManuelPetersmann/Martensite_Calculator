function [] = ave(solutions, i1, i2, prop_string, xi)

sol1 = solutions.array(i1);
sol2 = solutions.array(i2);

P1 = sol1.(prop_string);
P2 = sol2.(prop_string);

%[l1, l2, l3] = sorted_eig_vals_and_vecs( P1'*P1 )
%[l1, l2, l3] = sorted_eig_vals_and_vecs( P2'*P2 )

P_composite = xi* P1 + (1-xi)*P2;
(det(P_composite) - 1.023310863870360) * 100 % detB3

%% via polar decomposition
%[U1,R1] = polardecomposition(P1);
%[U2,R2] = polardecomposition(P2);

%F_c2 = rotation_average( R1, R2 )* ( xi* U1 + (1-xi)*U2 );
%det(F_c2)

%% without decomposing the shear into rotation and stretch and by using the rotation Q and ST

%F_c3 = rotation_average( sol1.Q, sol2.Q )* ( xi* sol1.ST + (1-xi)*sol2.ST );
%det(F_c3)


end