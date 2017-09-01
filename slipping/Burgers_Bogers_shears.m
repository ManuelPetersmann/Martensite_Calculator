function [ S ] = Burgers_Bogers_shears( ) %a_aust, a_mart )
%BURGERS_BOGERS_SHEARS 

a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 
% a_aust / a_mart

I = eye(3);
% The first shear is a_gamma/18 [1-2-1] (equals sqrt(6)/18*unit-vector) on every (11-1)_gamma plane 
% and the second shear is a_gamma/16 [-12-1] (equals sqrt(6)/16*unit-vector) on every (111)_gamma plane

a_aust/18. % = a_gamma/18
a_aust/16 % = a_gamma/16

%plane_families_P2 =     [ [1 1 0] ];
%direction_families_P2 = [ [1 1 1] ];

g1 = 1./18. ;
g2 = 1./16. ;
d1 = [1; -2; -1];
n1 = [1; 1; -1];
d2 = [1; -2; -1];
n2 = [1; 1; 1];

S1  = (d1  * n1') ;
S2  = (d2  * n2') ;

% (g1*S1) * (g2 * S2)  -- ich verstehe nicht warum aber dieser Term liefert
% immer genau 0-matrix
% S =  (I + g1 * S1) * (I + g2 * S2)
S =  (I + g1 * S1 + g2 * S2);
[eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(S,eye(3))
%eps_0

[U,R] = polardecomposition( S );
% Interstingly this deformation automatically has an eigenvalue of 1
eigs(U);

% define Bain-strain
eta1 = (a_mart/a_aust)*sqrt(2);
eta3 = a_mart / a_aust; % this is one form of three pcossible for the bain strain

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf 1.0 getuned
% werden. Die differenz ist also (n1-1).
B3 = [eta1 0    0   
       0  eta1  0
       0  0  eta3];
   
%S2 = inverse(B3)*S; % this gives the additional shear to B to obtain S_burgers_bogers   
%[eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(S2,eye(3));
%eps_0
% This turns out not to be an IPS ???!

end

