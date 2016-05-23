function [ fs ] = all_from_family_symmetries( f, SG )
% ALL_SYSTEMS_FROM_FAMILY 
% This function takes a family of directions e.g <0 0 1> or planes {0 0 1}
% as a Nx3 matrix and calculates all systems of the specified symmetry group SG

fs = [0 0 0];
k=0;
for i = 1:size(f,1)
    for j = 1:size(SG,3)
        k = k+1;
        % Drehung des Vektors bei festem Koordinatensystem
        fs(k,:) = SG(:,:,j) * f(i,:)'; % Drehung des Koordinatensystems R^-1*v
    end
end

end

