function [ strain ] = directional_strain( dir, F )
% given a colum vector dir and a deformation gradient
% this function calculates the strain in the given direction - see e.g. 
% Bhattacharya - Microstructures of Martensite p.24
strain = sqrt(dot(dir, (F'*(F*dir)))) -1. ;
end

