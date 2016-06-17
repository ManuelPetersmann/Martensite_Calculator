function [ theta_p ] = calc_misorientation_angle( plgamma, plalpha )
%CALC_MISORIENTATION_ANGLE Summary of this function goes here
%   Detailed explanation goes here

% calculation of angle between two vectors in 3D space 
% (e.g. normal vectors of two planes, in order to determine misorientation angle)
theta_p = acosd(dot(plgamma,plalpha)/(norm(plgamma)*norm(plalpha)));

end

