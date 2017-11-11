function [ R ] = axis_angle_to_rotmat( alpha, n )
% call as: rotmat_to_axis_angle( alpha, n )
% v_in - (row) vector to rotate
% alpha - angle to rotate in ° - degree
% n - around which should be rotated (right-handed system, positive sense)

if size(n'*n) ~= [3,3]
    %error('stupid matlab allows skalar + matrix - transpose vector')
    n = n';
end

alpha = deg2rad(alpha); % umrechnen von Grad aus Argument auf Radianten for cos und sin
n = n / norm(n);

% Kreuzproduktmatrix
nk = [ 0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0 ];
R = ( eye(3) *cos(alpha) +  ( 1. - cos(alpha))*(n'*n) + nk* sin(alpha) ); %   *  v_in'

% or equally writing v_in in each term respectively
%v_out2 = v_in*cos(alpha)  +  ( 1 - cos(alpha))* n*dot(n,v_in)  + cross( n , v_in ) * sin(alpha);
end
