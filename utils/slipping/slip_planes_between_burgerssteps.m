function [ stepwidth ] = slip_planes_between_burgerssteps( b, eps_s, plane_miller, lattice) 
% call: slip_planes_between_burgersstep( b, eps, plane_miller)
% b... Burgers vector of slip system (miller indizes)
% eps... shear magnitude of simple shear given below
% plane_miller... normal vector of slip system plane (miller indizes) 
% Given a simple shear of the form 
% S = eye(3) + eps  unit_slip_direction \otimes unit_slip_plane_normal
% calculates the average number of of slip_planes between Burgers steps "m"
% used in some formulations of crystallographic slip see e.g. Khachaturyans
% book: Theory of structural transfomrations in solids
% The function is based on the Intercept theorem (simple shear with unit
% vectors vs variable vectors) and the distance between lattice planes
% The fourth argument is a string identifying the Bravais lattice:
% Possible values are: 'cubic', 'tetragonal','hexagonal',rhombohedral',
% 'monoclinic','triclinic',

for i = 1:length(eps_s)
    % Note: in the case of the cubic lattice the lattice parameters fall out,
    % therefore they are not given here! TODO generalize burgersvector with
    % lattice parameters in this function (now first argument
    d = interplanar_distance( 1., plane_miller(i,1:3), lattice );
    stepwidth(i) = norm(b(i,1:3)) / (eps_s(i) * d);
    % Note:  slip density = (1/stepwidth)
end


end

