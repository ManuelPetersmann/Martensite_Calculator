function [strain] = write_strain_from_ST( filename, mat_array, large )
% call: write_strain_from_ST( file_name, mat_array, large )

if nargin < 3
    large = 0;
end

fileID = fopen(filename, 'a');
format_block_variant = '%7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \n'; %\t %7.4f \t %7.4f \t %7.4f \t \n';
for i=1: size(mat_array,3)
    % calculate strains (infinite deformation)
    if large==0
        strain(:,:,i) = 0.5*( mat_array(:,:,i)' + mat_array(:,:,i)) - eye(3);
    else
        strain(:,:,i) = 0.5*( mat_array(:,:,i)' * mat_array(:,:,i) - eye(3) );
    end
    fprintf( fileID,format_block_variant, get_order( strain(:,:,i) ) );
end
fclose('all');

    function six_entries = get_order( M )
        six_entries = [M(1,1) M(2,2) M(3,3) M(1,2) M(2,3) M(1,3)]; % Order definition as in my Zebfront Code % M(2,1) M(3,2) M(3,1)];
    end

end


