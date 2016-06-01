martensite = Martensite(); % creates martensite object
austenite = Base();

austenite.my_base =  [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
martensite.my_base = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
austenite.Bravais_type = 'cubic';
martensite.Bravais_type = 'cubic';
% austenite.Centering = 
% martensite.Centering = 'I';
a_aust = 3.6017264 % for 140Grad, 3.5975576 for 80Grad
a_mart = 2.8807346 % for 140Grad, 2.8790068 for 80Grad - check if something changes 
austenite.Lp = a_aust*[1 1 1];  % 3.5975576 %{1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
           % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
martensite.Lp = a_mart *[1 1 1];  % 2.8807346     

% Correspondance Matrix for B3 Bain
C_am = [0.5,  -0.5,  0.;
        0.5,   0.5,  0.;
        0.,    0.,   1.] ; 
%correspondance_matrix(austenite, martensite, m1, m2);

% define Bain-strain
n1 = (a_mart/a_aust)*sqrt(2);
n3 = a_mart / a_aust; % this is one form of three possible for the bain strain
% Der mittlere Eigenwert ist hier also n1. Dieser soll auf eins getuned
% werden. Die differenz ist also (n1-1).
B3 = [n1 0    0   
       0  n1  0
       0  0  n3];
martensite.B = B3;

det(B3)
%vars = martensite.variants()

cp = B3*C_am;

% assemble slip systems in alpha
[ ns, ds ] = independent_slipsystems_alpha();
% transform them to austenite


% highly symmetric mirror planes from bcc
% {001} family
m(1,:) = [0. 0. 1.];
m(2,:) = [0. 1. 0.];
m(3,:) = [1. 0. 0.];
% {011} family
m(4,:) = [0. 1. 1.];
m(5,:) = [1. 0. 1.];
m(6,:) = [-1. 0. 1.];
m(7,:) = [1. 1. 0.];
m(8,:) = [0. -1. 1.];
m(9,:) = [1. -1. 0.];

