clc
clear all;

a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 

Bain_and_Correspondence;

%% assemble slip systems in alpha
% since the shear is a substantial part of the transformation only 
% shear systems which are favorable in the b.c.c. lattices are considered. 
% the plane and direction families are {110}_alpha, {112}_alpha,
% <111>_alpha, <110>_alpha
plane_families_bcc =     [ [1 1 0]
                           [1 1 2] ];   % must be written with linebreak or ";" between vectors!                     
direction_families_bcc = [ [1 1 1]
                           [1 1 0] ];
                       
count_directions_extra = true;
                       
% find all possible combination (including different shear directions)
[martensite.slip_planes, martensite.slip_directions] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );
%[ ns_product, ds_product ] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );

plane_families_fcc =     [ [1 1 1] ];
direction_families_fcc = [ [1 1 0]; [1 1 2] ];
[austenite.slip_planes, austenite.slip_directions] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);
%[ ns_parent, ds_parent] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);

martensite.considered_plasticity = 3; % 1-mart, 2-aust, 3-both mart and aust slip systems
    
cpps_gamma = all_from_family_perms( [1 1 1] );
austenite.CPPs = cpps_gamma;

[ds, ns, S, slip_combinations] = shear_dyads(martensite, austenite, false); % assemble normed- shear_dyads

% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
% densest packed direction in austenite
% KS = u !!!!
us = all_from_family_perms( [1 1 0] ); %, false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
us = us / sqrt(2);

u = us(1,:)';

%%

%length(S)
%u = [1. 1. 1.]'

fid = fopen('equivalent_screws','a');
fprintf(fid,'%s %t %s %t %s',' s ',' b ',' m');
%for i = 1:length(u)
for i = 1:length(S)    
%    for j = i+1:length(S)
        
        [screw_syst] = equivalent_shear_screwdisloc( ds(i,1:3), ns(i,1:3) );
       
        fprintf(fid,'\n %s', [mat2str(ds(i,1:3)),'  ',screw_syst] );

        
%        u = us(i,:)';
%        u2 = B3 * ( eye(3) + (1./10)*S(:,:,j) ) * ( eye(3) + (1./20)*S(:,:,i) ) * u;
%        u2 = B3 * u;
%         
%         % calculate R such that u is unrotated
%         R = rotation_between_vectors( u2, u );
%            try
%            [angle, ax] = rotmat_to_axis_angle( R );
%            catch
%            end
%         
%         u2_mod = R * B3 * u;
%         
%         res = norm( u - u2_mod );
%         if  res < 0.02
%             i
%             j
%             res
%             ax
%         end
        
%    end
end

