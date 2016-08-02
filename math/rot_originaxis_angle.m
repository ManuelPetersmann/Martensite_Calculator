function Rot = rot_originaxis_angle( alpha, n )
% call: rot_originaxis_angle( alpha, n )
% rotation through origin axis, alternative to eulerangle rotation
% alpha is in degree!!! 

alpha = alpha*(pi/180); % grad to radians      
n1=n(1);
n2=n(2);
n3=n(3);
         
Rot = [(n1^2)*(1-cos(alpha))+cos(alpha)      n1*n2*(1-cos(alpha))-n3*sin(alpha)      n1*n3*(1-cos(alpha))+n2*sin(alpha)
        n2*n1*(1-cos(alpha))+n3*sin(alpha)   (n2^2)*(1-cos(alpha))+cos(alpha)         n2*n3*(1-cos(alpha))-n1*sin(alpha)
        n3*n1*(1-cos(alpha))-n2*sin(alpha)   n3*n2*(1-cos(alpha))+n1*sin(alpha)       (n3^2)*(1-cos(alpha))+cos(alpha)];
 
 end

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
