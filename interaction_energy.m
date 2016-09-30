function [ W_int ] = interaction_energy( epsilon1, epsilon2, n, E, nu )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lam = nu*E / ((1+nu)*(1-2*nu));
mu = E / (2*(1+nu));

n = n / norm(n);

% calculation of interaction matrix
% following C. Niclaeys et al., Int. J. Plast. 18 (2002), 1619-1647

% declaration and initialization of variables
F = zeros(3,3,3,3);
Q = zeros(3,3,3,3);
I = eye(3);

% calculation of interface operator F^{nm} between variants n and m,
% using the interface N^{nm} between the domains of variants n and m,
% see Eq.(12)
N = n * n';
%
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                F(i,j,k,l) = 0.5 * ( (I(i,k) - N(i,k)) * ( I(j,l) - N(j,l) ) ) + 0.5 * ( (I(j,k) - N(j,k)) * ( I(i,l) - N(i,l) ) ) ;
            end
        end
    end
end

% calculation of interface operator Q^{nm}, using F^{nm}
% used to relate the jump of internal stresses across interface
% between the phases n and m and the transformation strains across interf.,
% see Eq. (11)
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Q(i,j,k,l) = 2*mu * F(i,j,k,l) * lam / (lam+2*mu) * ( I(i,j) - N(i,j) ) * ( I(k,l) - N(k,l) );
            end
        end
    end
end

% calculation of a_n
a_n = 0.;
a_m = 0.;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                a_n = a_n + epsilon1(i,j)*epsilon1(k,l)*Q(i,j,k,l);
                a_m = a_m + epsilon2(i,j)*epsilon2(k,l)*Q(i,j,k,l);
            end
        end
    end
end
a_n = sqrt(a_n);
a_m = sqrt(a_m);

% calculation of interaction energy W_{int}, using Q^{ijkl} and a_n, a_m
% see Eq. (13)
W_int = 0.;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3              
                W_int = W_int + 0.5 / (a_n*a_m) * ( a_m * epsilon1(i,j) - a_n * epsilon2(i,j) ) * Q(i,j,k,l) * ( a_m * epsilon1(k,l) - a_n * epsilon2(k,l) );
            end
        end
    end
end


end
% calculation of interaction matrix H from W_{int}, relation Eq.(13)
