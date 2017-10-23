function write_lath_solutions(filename, permission, mart, reduced_solutions)
% write the (reduced,sorted) lath solutions

%FileID = fid
fid = fopen(filename,permission);
fmat = '[ %7.4f %7.4f %7.4f ]' ; 

if nargin < 4
    sols = mart.IPS_solutions;
else
    sols = reduced_solutions;
end

fprintf( fid,'\n %s \n\n',...
'###########################  LATH SOLUTIONS  ###########################');
fprintf( fid,'All variables are given in the coordinate system of the parent phase (austenite) all angles (theta_) are in degree! \n \n');
%
%% general information of IPS output data
% fprintf( fid,'Unphysical solutions with  delta_determinant_max > 0.1% are sorted out automatically')
fprintf( fid,'--- Invariant plane strain variables --- \n');
fprintf( fid,'%s \t\t ... \t %s \n','h ','habit plane normal vector');
fprintf( fid,'%s \t\t ... \t %s \n','d ','shape deformation direction');
fprintf( fid,'%s ... \t %s \n','eps_ips ','shape deformation magnitude');
fprintf( fid,'%s \t ... \t %s \n','ST ','full shape strain (including everything)');
fprintf( fid,'%s \t ... \t %s \n','LT ','lattice deformation (excluding "invisible" plastic deformation)');
fprintf( fid,'%s \t\t ... \t %s \n','Q ','rotation of inclusion to obtain an invariant plane strain');

fprintf( fid,'\n\n--- Slip system variables ---\n');
fprintf( fid,'%s \t\t ... \t %s \n','m ','slip plane normal');
fprintf( fid,'%s \t\t ... \t %s \n','s ','shear direction (Burgers vector in glide-plane)');
fprintf( fid,'%s \t ... \t %s \n','eps_s ','shear magnitude of slip system');
fprintf( fid,'%s... \t%s \n','stepwidth','averaged dislocation spacing in glide system (unit is: glide plane spacing).');
% TODO extend

keys = sols.selection_criteria.keys;
%values = sols.selection_criteria.values;
fprintf( fid,'\n\n--- Further selected active selection criteria variables ---\n');
for k = 1:length( keys )
    switch keys{k}
        case 'theta_CPPs' 
            fprintf( fid,'%s \t\t\t\t... \t %s \n','theta_CPPs',' maximum misorientation angle of CP relation');  
            % TODO generally between two plane families in parent and product phase specified by the user (e.g. close packed planes cpps) ')
        case 'theta_h_to_CPP' 
            fprintf( fid,'%s \t\t\t... \t %s \n','theta_h_to_cpp',' maximum misorientation angle between habit plane and nearest close packed plane {111}_gamma');
        case 'theta_KS_min'
            fprintf( fid,'%s \t\t\t...  %s \n','theta_KS_min',' maximum deviation angle of ideal Kurdjumov-Sachs direction parallelism ( <110>_aust || <111>_mart)');
        case 'theta_NW_min'
               fprintf( fid,'%s \t\t\t... %s \n','theta_NW_min',' maximum deviation angle of ideal Nishiyama-Wassermann direction parallelism ( <112>_aust || <110>_mart )');
        case 'delta_determinant_max'
            fprintf( fid,'%s \t... \t %s \n','delta_determinant_max','maximum deviation of theoretical volume change from Bain strain [%]');
        case 'theta_preferred_ILSdir_to_h'
            fprintf( fid,'%s ... %s \n','theta_preferred_ILSdir_to_h','maximum angle between preferred invariant direction (invariant line strain) and habit plane (0° if vector is in habit plane)');
    end
end


%% function for writing slip system information
    function [s,m] = write_slipsys_strings( sol )
        for j = 1:size(sol.slip_normal_plane_vec,1)
            m = sprintf('(%d, %d, %d)_%s', sol.slip_normal_plane_vec(j,1:3), phase(sol.slip_normal_plane_vec(j,4) ) );
            s = sprintf('[%d, %d, %d]_%s', sol.shear_direction(1,1:3),       phase(sol.shear_direction(j,4) ) );
            fprintf( fid, ' m = %s \t s = %s \t eps_s = %5.4f \n', m, s, sol.eps_s(j) );
        end
        %
        function phase_string = phase(identifier)
            if identifier == 88
                phase_string = 'aust';
            end
            if identifier == 99
                phase_string = 'mart';
            end
        end
    end

%% loop over solutions for writing
fprintf( fid,'\n\n ---Solutions--- ');
%
% solutions sorted after 
if strcmp( mart.IPS_solutions.sorted_after,'unsorted')
    fprintf(fid,'\n %s', 'Solutions are unsorted');
else
    fprintf(fid,'\n %s', ['Solutions have been sorted ascendingly for property: ', mart.IPS_solutions.sorted_after] );
end
% length( sols.array )
for i=1:length( sols.array )
    
    s = sols.array(i);
    
    fprintf( fid,'\n\n\n %s %d \n','solution ID = ', s.id);
    fprintf( fid,[' h = ',fmat,'\t d = ',fmat,'\t  eps_ips = %5.4f \n'], s.h, s.d,s.eps_ips);
    write_slipsys_strings( s );
    %fprintf( fid,'%s \t \t \t %s \t \t \t %s','ST =','LT =','Q =');
    fprintf( fid,'%s \n','F1 = R*B*S2*S1 = ');
    fprintf( fid,[fmat,'\n'],s.F1);
    
    if ~isempty( s.added_props )
        sk = s.added_props.keys;
        % write info on selection criteria if any are specified
        for k = 1:length( sk ) % or keys
            switch sk{k}
                case 'theta_CPPs' % (minimum) misorientation angle of between two close packed planes (cpps) in parent and product phase
                    fprintf( fid,' %s %5.4f%s (%d, %d, %d)_%s \n','theta_CPPs   = ',s.added_props('theta_CPPs'),'° to ',s.added_props('closest_CPPs'),'aust ');
                    %
                case 'theta_h_to_CPP' % (minimum) misorientation angle between habit plane and nearest close packed plane
                    fprintf( fid,' %s %5.4f%s (%d, %d, %d)_%s \n','theta_h_to_CPP = ',s.added_props('theta_h_to_CPP'),'° to ',s.added_props('closest_h_to_CPP'),'aust ');
                    %
                case 'theta_KS_min'
                    fprintf( fid,' %s %5.4f%s [%d, %d, %d]_%s \n','theta_KS_min = ',s.added_props('theta_KS_min'),'° to ',s.added_props('closest_KS'),'aust ');
                    %
                case 'theta_NW_min'
                    fprintf( fid,' %s %5.4f%s [%d, %d, %d]_%s \n','theta_NW_min = ',s.added_props('theta_NW_min'),'° to ',s.added_props('closest_NW'),'aust ');
                    %
                    %             case 'delta_determinant_max'
                    %                 fprintf( fid,' %s = %5.4f \n','delta_determinant_max = ', s.delta_determinant_max);
                    %
                case 'theta_preferred_ILSdir_to_h'
                    fprintf( fid,' %s %5.4f%s [%d, %d, %d]_%s \n','theta_preferred_ILSdir_to_h = ',s.added_props('theta_preferred_ILSdir_to_h'),'° to ',s.added_props('closest_ILSdir_to_h'),'_aust ');
            end
        end
    end
end


fclose(fid);



end


