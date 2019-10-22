function [ UVJW ] = cubic_to_hexagonal_indices( uvw )
% takes a vektor of cubic miller indices and calculates the corresponding 
% hexagonal miller indizes

u = uvw(1);
v = uvw(2);
w = uvw(3);

U = (2*u-v)/3;
V = (2*v- u)/3;
J = -(U + V);
W = w;

UVJW = [U, V, J, W];
end

