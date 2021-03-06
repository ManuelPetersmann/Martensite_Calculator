
%% Set up coordinate systems and lattice parameters
martensite = Martensite(); % creates martensite object
austenite = Base();

% conventional bases
austenite.my_base =  [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
martensite.my_base = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
austenite.Bravais_type  = 'cubic';
martensite.Bravais_type = 'cubic';
% austenite.Centering = 
martensite.Centering = 'I';
a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 
austenite.Lp = a_aust*[1 1 1];  % 3.5975576 % {1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
                                              % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
martensite.Lp = a_mart *[1 1 1];  % 2.8807346   

%% Correspondence Matrix for B3 Bain
% e_mart = C_am * e_aust; infinetively many choices. Take one for slip
% system transformation....

martensite.C_am = [0.5,  -0.5,  0.;
                   0.5,   0.5,  0.;
                   0.,    0.,   1.]; 
%correspondance_matrix(austenite, martensite, m1, m2);

% define Bain-strain
eta1 = (a_mart/a_aust)*sqrt(2);
eta3 = a_mart / a_aust; % this is one form of three pcossible for the bain strain

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf 1.0 getuned
% werden. Die differenz ist also (n1-1).
B3 = [eta1 0    0   
       0  eta1  0
       0  0  eta3];
martensite.U = B3;

display(['Volume change in percent is:', num2str( ( det(B3)-1)*100 ) ] );
%vars = martensite.variants()

% moved into martensite class
%cp = B3*C_am;