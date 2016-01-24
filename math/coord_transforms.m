% This file provides several functions for COORDINATE TRANSFORMATIONS 
% such as rotating (passive rotation), Mirroring, scaling, shearing

function funcs = coord_transforms()
funcs.eulerAngleRotation =   @eulerAngleRotation;
funcs.rot_mat_axis =         @rot_mat_axis;
funcs.rot_vec_axis_angle =   @rot_vec_axis_angle;
funcs.rot_vec_axis_angle2 =  @rot_vec_axis_angle2;
end

% This class provides several functions for rotations
% of vectors and tensors
% (active and passive)


% TODO specifiy other conventions, i.e. order of rotations and different
% axis using euler angles

function [ v_i_out ] = eulerAngleRotation( phi1, theta, phi2, v_i_in)
% This function takes the three euler Angles after Bunges definition and
% first computes the rotation matrix with the z, x', z'' (= z' since the first rotation leaves z invariant)
% the function takes as input the three euler Angles ( in order of
% rotation) phi1, theta, phi2 and a 3x3 matrix that has the basis vectors
% of the cartesian coordinate system that should be transformed as colums.
% The output of the function is a matrix that has the rotated basis vectors
% as colums

Rz = [cos(phi2) sin(phi2) 0;
     -sin(phi2) cos(phi2) 0;
         0         0      1];

Rx_strich = [ 1      0          0;
              0  cos(theta) sin(theta);
              0 -sin(theta) cos(theta)];
          
Rz_2strich = [ cos(phi1)  sin(phi1)  0;
              -sin(phi1)  cos(phi1)  0;
                    0         0      1];
                
Rot = Rz * Rx_strich * Rz_2strich;

v_i_out = zeros(3);

for i = 1:3
    % note that the rotation matrix has to be (only {because it is
    % orthogonal})transposed here!
    v_i_out(:,i) = Rot' * v_i_in(:,i);
end

end

%-------------------------------------------------------------------

function[ rotmat ] = rot_mat_axis( n, alpha )
% call as: rot_vec_axis_angle(v_in, alpha, n )
% v_in - (row) vector to rotate
% alpha - angle to rotate in Degree °
% n - vector through the origin around which should be rotated
%(right-handed system, positive sense)
alpha = deg2rad(alpha); %convert Grad to radians
n1=n(1) / norm(n,2);
n2=n(2) / norm(n,2);
n3=n(3) / norm(n,2);

rotmat = [(n1^2)*(1-cos(alpha))+cos(alpha)      n1*n2*(1-cos(alpha))-n3*sin(alpha)      n1*n3*(1-cos(alpha))+n2*sin(alpha)
    n2*n1*(1-cos(alpha))+n3*sin(alpha)   (n2^2)*(1-cos(alpha))+cos(alpha)         n2*n3*(1-cos(alpha))-n1*sin(alpha)
    n3*n1*(1-cos(alpha))-n2*sin(alpha)   n3*n2*(1-cos(alpha))+n1*sin(alpha)       (n3^2)*(1-cos(alpha))+cos(alpha)];
end

%-------------------------------------------------------------------

function[ v_out ] = rot_vec_axis_angle( v_in, alpha, n )
rot_mat_axis( alpha, n )
v_out = Rot * v_in';
end

%-------------------------------------------------------------------

function [ v_out ] = rot_vec_axis_angle2( v_in, alpha, n )
% call as: rot_vec_axis_angle(v_in, alpha, n )
% v_in - (row) vector to rotate
% alpha - angle to rotate
% n - around which should be rotated (right-handed system, positive sense)
alpha = deg2rad(alpha); 
n = n / norm(n);
% crossproductmatrix
nk = [ 0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0 ];
v_out = ( ( 1 - cos(alpha))*(n'*n) + eye(3)*cos(alpha) + nk*sin(alpha) ) * v_in';
end