function [E] = El_from_F( F )
% calculates green lagrange / large strain E from Deformation gradient F

E = 0.5*(F'*F - eye(3) );
end