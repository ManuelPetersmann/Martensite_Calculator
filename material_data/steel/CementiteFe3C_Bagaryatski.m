cementite = Base([1. 0. 0.; 0. 1. 0.; 0. 0. 1.]);
ferrite = Base([1. 0. 0.; 0. 1. 0.; 0. 0. 1.;]);
cementite.Bravais_type = 'orthorhombic';
ferrite.Bravais_type = 'cubic';
%ferrite.Point_group
%cementite.Centering = 
ferrite.Centering = 'I';
ferrite.Lp(1:3) = 2.8662; %{1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
           % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
cementite.Lp(1:3) = [4.5241, 5.0883, 6.7416];           

% Bagaryatski OR
% directions in cementite
m1 = [1 0 0;
      0 1 0;
      0 0 1];
% directions in ferrite
m2 = [0 1  2;
     -1 -1 1;
      1 -1 1];

or = lattice_correspondance();
cp = or.correspondance_matrix(cementite, ferrite, m1, m2);
