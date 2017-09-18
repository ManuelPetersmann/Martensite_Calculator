function ave_rot = rotation_average( R1, R2 )
% call: rotation_average( R1, R2 )
% where R1 and R2 are rotation matrices
% Averaging of two rotation matrices due to Moahker 2002
% Also see Cayron 2016 Crystallography 557... 
% This way e.g. Distortion matrices of Martensite with the same Bain axis 
% (same Bain strain in the Polar decomposition) can be averaged by
% averaging their rotation matrices such that the determinant of the
% deformation is conserved. 
% Note that the determinant is not a linear function of matrices

R1 * sqrtm(R1' * R2);
ave_rot = R2 * sqrtm(R2' * R1);


end

