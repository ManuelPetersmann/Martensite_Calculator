clc

a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 

Bain_and_Correspondence;
    
cpps_gamma = all_from_family_perms( [1 1 1] );
austenite.CPPs = cpps_gamma;

% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
% densest packed direction in austenite
% KS = u !!!!
u = all_from_family_perms( [1 1 0] ); %, false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
u = u / sqrt(2);

for i = 1:length(u)
    % calculate R such that u is invariant line
    u(i,:)'
    u2 = B3 * u(i,:)'
    %
    ax = cross( u2, u(i,:) ); % sense of rotation is from KS(1,:)' to u2
    ax = ax / norm(ax);
    %
    ang = acosd( dot ( u(i,:), u2 ) / ( norm(u(i,:)) * norm(u2) ) )
    %signed_angle_from_rotmatrix( R )
    if abs(ang) > 1.e-4
        % no unrotated candidate for ILS!
        continue
    end
    %
    R = axis_angle_to_rotmat( ang, ax ); % n = all KS!
    % either I change the direction of rotation or I take the negative of each
    % vector! (the latter is used now!)
    
    %B3 * R * KS(1,:)'
    KS_unrot(i,:) = R * B3 * u(i,:)';
    abs(1. - norm( KS_unrot(i,:) ) )
    
end

%%
% [ y1, y2, y3, e1, e2, e3] = sorted_eig_vals_and_vecs( Ct )
% % This function takes a matrix and returns its eigenvalues in ascending 
% % order as well as the corresponding (normalized!) eigenvectors 
% 
% [V,D] = eig( Ct );


%while
    
    
    
%end