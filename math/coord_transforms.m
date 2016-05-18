% This file provides several functions for COORDINATE TRANSFORMATIONS 
% such as rotating (passive rotation), Mirroring, scaling, shearing
% If not stated otherwise the variable "input" is a 1x3 colon vector, 3x3
% matrix or 3x3x3x3 tensors and the resulting variable "out" has the same
% sizes

% The principle of frame indifference is that the free energy does not
% change for a passive rotation in the symmetry group

function funcs = coord_transforms()
funcs.rot =                  @rot;
funcs.eulerAngleRotation =   @eulerAngleRotation;
funcs.rot_mat_axis_angle =   @rot_mat_axis_angle;
funcs.rot_axis_angle =       @rot_axis_angle;
%funcs.rot_vec_axis_angle2 =  @rot_vec_axis_angle2;
end

%-------------------------------------------------------------------
function out = rot(input, rot_mat)
if size(input) == [1, 3] %input must be a colon vector
    % Note that generally instead of the transpose here is an inverse, but
    % since for a rotatoin R*R'=I holds its the same
    out = rot_mat' * input; % passive rotation
elseif size(input) == [3, 3]
    out = zeros(3);
    for i = 1:3
        out(:,i) = rot_mat' * input(:,i);
    end
elseif size(input) == [3, 3, 3, 3]
    % TODO implement for 4th order tensor
end
end
%-------------------------------------------------------------------

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
rot_mat = rot(input, rot_mat);
end
%-------------------------------------------------------------------
function rot_mats = rot_mat_axis_angle( n, angles )
% Rotation matrices around an axis given by a normalized vector n through the origin
% given a vector or a list of vectors (3xN) n and an angle or a list of angles (1xN) 
% returns a 3x3 or 3x3xN array of rotation matrices associated to them
% call as: rot_mat_axis_angle( n , alpha )
% n (row)-vector(s) to rotate, alpha - angle(s) to rotate (in Degree °)
%
if size(angles) == [1,1]
    angles = angles*ones(size(n, 2));
end
if size(n,2) ~= size(angles,2)
    error('vectorlist and angles must be the same size')
end
%
rot_mats = zeros(3,3,size(n, 2));
for i=1:size(n,2)
    alpha = deg2rad(angles(i)); %convert Grad to radians
    n1 = n(1,i) / norm(n(:,i),2);
    n2 = n(2,i) / norm(n(:,i),2);
    n3 = n(3,i) / norm(n(:,i),2);
    %
    rot_mats(:,:,i) = [(n1^2)*(1-cos(alpha))+cos(alpha)      n1*n2*(1-cos(alpha))-n3*sin(alpha)      n1*n3*(1-cos(alpha))+n2*sin(alpha)
        n2*n1*(1-cos(alpha))+n3*sin(alpha)   (n2^2)*(1-cos(alpha))+cos(alpha)         n2*n3*(1-cos(alpha))-n1*sin(alpha)
        n3*n1*(1-cos(alpha))-n2*sin(alpha)   n3*n2*(1-cos(alpha))+n1*sin(alpha)       (n3^2)*(1-cos(alpha))+cos(alpha)];
end
end
%-------------------------------------------------------------------
function out = rot_axis_angle( input, alpha, n )
rot_mat_axis( alpha, n );
out = rot(input, rot_mat);
end
%-------------------------------------------------------------------

% function [ v_out ] = rot_vec_axis_angle2( v_in, alpha, n )
% % call as: rot_vec_axis_angle(v_in, alpha, n )
% % v_in - (row) vector to rotate
% % alpha - angle to rotate
% % n - around which should be rotated (right-handed system, positive sense)
% alpha = deg2rad(alpha); 
% n = n / norm(n);
% % crossproductmatrix
% nk = [ 0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0 ];
% v_out = ( ( 1 - cos(alpha))*(n'*n) + eye(3)*cos(alpha) + nk*sin(alpha) ) * v_in';
% end