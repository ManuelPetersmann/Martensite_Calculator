function [ eps ] = normed_shearmag_from_planes_between_burgerssteps( b, stepwidth, plane_miller ) %, lattice) - up do now d only for cubic lattice!
% call: slip_planes_between_burgersstep( b, eps, plane_miller)
% b... Burgers vector of slip system (miller indizes)
% stepwidth ... average number of of slip_planes between steps in interface
% or equally averge -"- between dislocations in the bulk
% used in some formulations of crystallographic slip see e.g. Khachaturyans
% book: Theory of structural transfomrations in solids
% plane_miller... normal vector of slip system plane (miller indizes)
% Given a simple shear of the form 
% S = eye(3) + eps  unit_slip_direction \otimes unit_slip_plane_normal
% calculates its shear magnitude eps given a.m. input.
% The function is based on the Intercept theorem (simple shear with unit
% vectors vs variable vectors) and the distance between lattice planes

% I DO NOT NEED THIS - IF I REPLACE M WITH EPS IT is valid as well since
% m propto eps - hence the function slip_planes_between_burgersvec is
% sufficient for the transformation in both "directions"

for i = 1:length(stepwidth)
    % TODO add other crystal systems! and generalize Burgers vector with
    % lattice parameters in this function
    d(i) = 1/norm(plane_miller(i,:));
        
    eps(i) = norm(b(i,:)) / (stepwidth*d);
    
    % m = norm(b) / (eps * d);
end

end

