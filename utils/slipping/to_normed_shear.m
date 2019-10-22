function [ eps ] = to_normed_shear( b_miller, stepwidth, plane_miller, lattice)
% call: slip_planes_between_burgersstep( b, eps, plane_miller)
% b... Burgers vector of slip system: miller indizes in austenite or cp*b_miller in martensite
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
% Formally this function was called ...planes_between_burgerssteps (but
% this was based on the non-general assumption that the shear is normal to
% the inital interface

% I DO NOT NEED THIS - IF I REPLACE M WITH EPS IT is valid as well since
% m propto eps - hence the function slip_planes_between_burgersvec is
% sufficient for the transformation in both "directions"

if nargin < 4
    lattice = 'cubic';
end

%% TODO - generalize to other lattices - now only valid for cubic ones!
    % PET: correction 21.11.17 - corrected here factor 0.5 due to
    % Burgesvector of form a/2 * norm(b) = sqrt(3)*a /2 or sqrt(2)*a / 2
    % PET: 2.3.18: removed factor 1/2 again since it is only for one
    % specific shear direction - what must be used are the miller indices
    % of fcc or C_am * B * millers_bcc = cp * millers_bcc 

for i = 1:length(stepwidth)
    % TODO add other crystal systems! and generalize Burgers vector with
    % lattice parameters in this function - instead of one then below there
    % must be an array of lattice parameters!!!
    d(i) = interplanar_distance( 1., plane_miller(i,1:3), lattice ); %1/norm(plane_miller(i,:));
        
    eps(i) = norm(b_miller(i,:)) / (stepwidth*d);
    
    % m = norm(b) / (eps * d);
end

end

