function write_input_parameters(filename,permission,mart, aust)
%format long
%FileID = fid
fid = fopen(filename,permission);

% specify print formats
fmat = '[ %7.4f %7.4f %7.4f ]\n' ; %    [ %7.4f %7.4f %7.4f ]\n    [ %7.4f %7.4f %7.4f ]\n \n';
fplanefam = '{%d, %d, %d} \t\n';
fdirfam   = '<%d, %d, %d> \t\n';

% write program name
fprintf(fid, '%s \n','f.c.c - b.c.c highly dislocated, hierarchical - MARTENSITE-CALCULATOR - lath / block invariant plane strain (IPS) microstructure statistics');

% write date
fprintf(fid, '%s \n', ['file created: ',datestr(datetime('now'))] );
fprintf(fid, '\n %s \n',...
'##########################  INPUT PARAMETERS  ##########################');

% write lattice parameters
fprintf(fid, '\n %s \n','lattice parameters:');
fprintf(fid,' austenite  = %7.4f \t ',aust.Lp(1)  );
fprintf(fid,' martensite = %7.4f \n ',mart.Lp(1) );

% specific structural stretch tensor
fprintf(fid, '\n %s \n','specific Bain:');
fprintf(fid, fmat, mart.U);

% write correspondence matrix (coordinate trans)
fprintf(fid, '\n %s \n','Correspondence matrix (coordinate transformation for specific Bain):');
fprintf(fid, fmat, mart.my_base);

% write Bases
fprintf(fid, '\n %s \n','austenite basis:');
fprintf(fid, fmat, mart.my_base);
fprintf(fid, '\n %s \n','martensite basis');
fprintf(fid, fmat, aust.my_base);

% write slip plane/direction families of phases
if (mart.considered_plasticity == 1 || mart.considered_plasticity == 3)
    fprintf(fid, '\n %s \n','crystallographic families in martensite to build slip systems');
    fprintf(fid, fplanefam, mart.slip_plane_families );
    fprintf(fid, fdirfam,   mart.slip_dir_families);
end
if (mart.considered_plasticity == 2 || mart.considered_plasticity == 3)
    fprintf(fid, '\n %s \n','crystallographic families in austenite to build slip systems');
    fprintf(fid, fplanefam, aust.slip_plane_families);
    fprintf(fid, fdirfam,   aust.slip_dir_families);
end
fprintf(fid, '\n %s \n \t %d \n','resulting possible combinations of double shears (both slip directions, nchoosek, k=2): ', mart.IPS_solutions.slip_combinations );

fclose(fid);
end




        
