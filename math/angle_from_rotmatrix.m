function [ angle ] = angle_from_rotmatrix( R )
%ANGLE_FROM_ROTMATRIX Summary of this function goes here
%   Detailed explanation goes here

angle = acos( (trace(R) - eye(3)) / 2. );

end

