function write_opt_block_solutions(filename, permission,  block_sols)
% write the (reduced,sorted) lath solutions

%FileID = fid
fid = fopen(filename,permission);
fmat = '[ %7.4f %7.4f %7.4f ]' ; 

fprintf( fid,'\n %s \n\n',...
'###############  CALCULATION OPTIONS FOR OPTIMIZED BLOCKS  ###############');

% calculation method
fprintf(fid,'%s \n \t %s \n\n','Calculation method: ', block_sols.calculation_method ); 

fprintf( fid,' --- Active criteria allowing for mixing of lath solutions---\n');
bsmt = block_sols.mixing_tolerances.keys;
for k = 1:length( bsmt )
    switch bsmt{k}
        case 'theta_intersec_cpdir' % (minimum) misorientation angle of between two close packed planes (cpps) in parent and product phase
            fprintf(fid,'%s \n', 'Maximum allowed angle between intersection line of habit planes (largest dimension of laths) and close packed direction' );
            fprintf( fid,' %s %5.4f \n','theta_intersec_cpdir   = ',block_sols.mixing_tolerances('theta_intersec_cpdir') );
            %
        case 'theta_hps' % (minimum) misorientation angle between habit plane and nearest close packed plane
            fprintf(fid,'%s \n','Maximum allowed angle between habit planes (second largest dimension of laths)');
            fprintf( fid,' %s %5.4f \n','theta_hps = ',block_sols.mixing_tolerances('theta_hps') );
    end
end

% COULD BE EXTENDED select for - prop_to_optimize - TO OPTIMIZE OTHER NORMS, 
% e.g. Frob_green_lagrange - WHEN THERE IS NO VOLUME CHANGE... 

fprintf( fid,'\n\n %s \n\n',...
'#######################  OPTIMIZED BLOCK SOLUTIONS  #######################');

%% general information of Block output data
fprintf( fid,'%s \n','F_comp ... x * F1_lath1 + (1-x) * F1_lath2   ');
fprintf( fid,'%s \n','_composite_block ... Invariant plane strain properties build with Fcomp.');
%fprintf( fid,'%s \n','
fprintf( fid,'%s \n','x_eps ... optimized phase fractions [x, 1-x] of minimum euklidean norm of d_comp: linear mixture of shape strain vectors "d" of lath solutions.');
fprintf( fid,'%s \n','x_dis ... optimized phase fractions [x, 1-x] of minimum frobenius norm of:  F_comp - I              ( displacement gradient )');
fprintf( fid,'%s \n','x_gl  ... optimized phase fractions [x, 1-x] of minimum frobenius norm of:  F_comp^T * F_comp - I   ( Green-Lagrange tensor )');
%fprintf( fid,'---Invariant plane strain variables--- \n');
%fprintf( fid,'%s \t\t ... \t %s \n','h ','habit plane normal vector');

% solutions sorted after 
% if strcmp( mart.IPS_solutions.sorted_after,'unsorted')
%     fprintf(fid,'\n %s', 'Solutions are unsorted');
% else
%     fprintf(fid,'\n %s', ['Solutions have been sorted ascendingly for property: ', mart.IPS_solutions.sorted_after] );
% end

for i=1:length( block_sols.array )
    
    s = block_sols.array(i);
    
    fprintf( fid,'\n\n\n %s [%d , %d]','used lath solutions (IDs) = ', s.lath_id_pair);
    %
    if ~isempty( s.tolerances )
        sk = s.tolerances.keys;
        % write info on selection criteria if any are specified
        for k = 1:length( sk ) % or keys
            switch sk{k}
                case 'theta_intersec_cpdir' % (minimum) misorientation angle of between two close packed planes (cpps) in parent and product phase
                    fprintf( fid,'\t %s %5.4f ','theta_intersec_cpdir   = ',s.tolerances('theta_intersec_cpdir') );
                    %
                case 'theta_hps' % (minimum) misorientation angle between habit plane and nearest close packed plane
                    fprintf( fid,'\t %s %5.4f ','theta_hps = ',s.tolerances('theta_hps') );
            end
        end
    end
    %
    fprintf( fid,'\n');
    fprintf( fid,[' h_composite_block = ',fmat,'\t d_composite_block = ',fmat,'\n'], s.h, s.d);
    fprintf( fid,' eps_ips_composite_block = %5.4f \n',s.eps_ips);
    %
    fprintf( fid,' x_eps = [%3.4f , %3.4f] \t eps_net = %5.4f \n', s.x_eps(1), s.x_eps(2), s.eps_net);
    %fprintf( fid,'%s \n','F_comp_eps = ');
    %fprintf( fid,[fmat,'\n'],s.F_comp_eps);
    %
    s.x_dis
    s.x_gl
    %
    fprintf( fid,' x_dis = [%3.4f , %3.4f] \t frob_opt_displacement_grad = %5.4f \n',  s.x_dis, s.frob_opt_displacement_grad);
    s.x_gl
    fprintf( fid,' x_gl  = [%3.4f , %3.4f] \t frob_opt_green_lagrange =    %5.4f \n',  s.x_gl,  s.frob_opt_green_lagrange);    
end


% block_solutions.array(1).lath_solution_pair(1)
% block_solutions.array(1).lath_solution_pair(2)

fclose(fid);



end


