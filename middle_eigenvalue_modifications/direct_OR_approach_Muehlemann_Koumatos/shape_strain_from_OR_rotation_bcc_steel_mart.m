function [ ST ] = shape_strain_from_OR_rotation_bcc_steel_mart( r , omega, eta1, eta3  )
% given the average orientation relation in form of a rotation around the
% axis r with an angle omega and the entries eta1 and eta3 of the Bain strain,
% computes the according transformation shape strains. The function assumes 
% an -45 degree rotation around the Bain axis to relate the fcc and bcc systems
% (Correspondence matrix CP). - See Mühlemann 2016

OR = rot_originaxis_angle( omega , r );

% identify bain axis ba
ba_idx = 0.94;
for i = 1 : 3
    if r(i) > ba_idx
        ba_idx = i;
    end 
end

if ba_idx == 0.94
    error('high deviation from Bain axis are you sure this is a steel OR?')
end

ba = [0., 0., 0.];
ba(ba_idx) = 1. ;

% define correspondence matrix
CP = rot_originaxis_angle( -45. , ba );

Ri = (inverse(CP) * OR)';

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf 1.0 getuned
% werden. Die differenz ist also (n1-1).
B = [eta1 0    0   
       0  eta1  0
       0  0  eta1];
B(ba_idx,ba_idx) = eta3;   
   
ST = Ri * B;

end