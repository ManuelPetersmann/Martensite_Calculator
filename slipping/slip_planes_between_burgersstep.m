function [ m ] = slip_planes_between_burgersstep( b, eps, plane_miller ) %, lattice) - up do now d only for cubic lattice!
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

% TODO add other crystal systems!
d = 1/norm(plane_miller);

m = norm(b) / (eps * d);


end

