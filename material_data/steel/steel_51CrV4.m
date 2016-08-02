martensite = Martensite(); % creates martensite object
austenite = Base();


austenite.my_base =  [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
martensite.my_base = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
austenite.Bravais_type = 'cubic';
martensite.Bravais_type = 'cubic';
% austenite.Centering = 
% martensite.Centering = 'I';
%%% lattice constants in nm
%%% --> "Modeling the role of external stresses on the 
%%%      austenite-to-bainite phase transformation in 51CrV4 steel",
%%%      M.C. Uslu, D. Canadinc, H.-G. Lambers, S. Tschumak, H.J. Maier,
%%%      Modelling Simul. Mater. Sci. Eng. 19 (2011) 045007, (17pp)
%%%      for 51CrV4: a_gamma = 0.3605   und   a_alpha = 0.2884 in nm
% austenite:  face-centered-cubic (fcc)
a_aust = 3.605 % in Angström
% martensite: body-centered-cubic (bcc)
a_mart = 2.884 % in Angström


austenite.Lp = a_aust*[1 1 1];  
martensite.Lp = a_mart *[1 1 1];

% Correspondance Matrix for B3 Bain
C_am = [0.5,  -0.5,  0.;
        0.5,   0.5,  0.;
        0.,    0.,   1.]; 
    
det(B3)    
% % Correspondance Matrix for B3 Bain (other basis)
% C_am = [ 0.5,   0.5,  0.;
%         -0.5,   0.5,  0.;
%          0.,    0.,   1.];       
%correspondance_matrix(austenite, martensite, m1, m2);

% define Bain-strain
n1 = (a_mart/a_aust)*sqrt(2);
n3 = a_mart / a_aust; % this is one form of three possible for the bain strain

% temprary: rounding of values for direct comparison with results from paper
% n1 = round ( ((a_mart/a_aust)*sqrt(2))*1000 )/1000;
% n3 = round ( (a_mart / a_aust)*1000 )/1000;

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf eins getuned
% werden. Die differenz ist also (n1-1).
B3 = [n1 0    0   
       0  n1  0
       0  0  n3];
martensite.B = B3;

det(B3);
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