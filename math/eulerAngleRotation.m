% TODO specifiy other conventions, i.e. order of rotations and different
% axis using euler angles
function rot_mat = eulerAngleRotation( phi1, theta, phi2, input )
% This function takes the three euler Angles after Bunges definition and
% first computes the rotation matrix with the z, x', z'' (= z' since the first rotation leaves z invariant)
% the function takes as input the three euler Angles ( in order of
% rotation) phi1, theta, phi2 and a 3x1 vector, 3x3 matrix or 3x3x3x3 tensor 
% that has the basis vectors of the cartesian coordinate system that should
% be transformed as colums. The output of the function is a matrix that has
% the rotated basis vectors as colums
%
Rz = [cos(phi2) sin(phi2) 0;
     -sin(phi2) cos(phi2) 0;
         0         0      1];
%
Rx_strich = [ 1      0          0;
              0  cos(theta) sin(theta);
              0 -sin(theta) cos(theta)];
%          
Rz_2strich = [ cos(phi1)  sin(phi1)  0;
              -sin(phi1)  cos(phi1)  0;
                    0         0      1];
%
rot_mat = Rz * Rx_strich * Rz_2strich;
%
%rot_mat = rot(input, rot_mat);
end
