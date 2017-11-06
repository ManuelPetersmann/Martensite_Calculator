parent = Martensite();
product = Base();
parent.Bravais_type = 'cubic';
product.Bravais_type = 'tetragonal';
parent.Centering = 'P';
product.Centering = 'I'; %body centered tetragonal
parent.Lp(1:3) = 2.869; 
           % {1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
           % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
product.Lp(1:3) = 3.56;
n3 = 0.806;
n1 = 1.140;

parent.Point_group;           
parent.E;
 
% KS OR
% {100}_aust || {100}_mart
% <100>_aust || <110>_mart

% Khattachuryan, Qi:
% (111)_aust || (011)_mart
% <1-10>_aust || <11-1>_mart
% Cayron:
% (111)_aust || (110)_mart
% [-111]_aust || [100]_mart

C = [0.5 -0.5 0; % should be 1/sqrt(2)
     0.5  0.5 0;
     0.   0.  1]
 
% F = C * B

parent.B = []