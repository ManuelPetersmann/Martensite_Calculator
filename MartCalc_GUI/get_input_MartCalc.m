% Script for reading input from MartCalc-GUI

% lattice parameter for fcc lattice
% example: MARVALX12 - a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
if(num2str(handles.lc_edtxt_aust_val.String) > 0 )
    a_aust = str2num(handles.lc_edtxt_aust_val.String)
else
    error('No reasonable input for fcc lattice parameter!');
end


% lattice parameter for bcc lattice
% example: MARVALX12 - a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 
if(num2str(handles.lc_edtxt_mart_val.String) > 0 )
    a_mart = str2num(handles.lc_edtxt_mart_val.String)
else
    error('No reasonable input for bcc lattice parameter!');
end
 


% criterion 1: Minimum slip plane density
if(handles.asc_status(1) > 0)
    g_min = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(1)).Children(2).String);
end
      
% Criterion 2: Maximum shape strain
if(handles.asc_status(2) > 0)
    eps_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(2)).Children(2).String);
end

% Criterion 3: Maximum misorientation of CPPs {110}_alpha and {111}_gamma
if(handles.asc_status(3) > 0)
    theta_p_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(3)).Children(2).String);
end

% Criterion 4: Maximum misorientation of block HP to {111}_gamma
% note 557 is 9.4° from 111 ! therefore this high tolerance!
% Angle between 111 and 557 habit plane
% acos( dot([1. 1. 1.], [5. 5. 7.])/(sqrt(3)*sqrt(99) ) ) = 9.4 degree
if(handles.asc_status(4) > 0)
    theta_n_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(4)).Children(2).String);
end

% Criterion 5: Maximum deviation of determinant det(F) of transformation
if(handles.asc_status(5) > 0)
    delta_determinant_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(5)).Children(2).String);
end

% Criterion 6: Maximum deviation from KS OR
% The peak of OR distribution is normally between KS and NW and
% these two are 5.25 apart - hence these tolerances
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
if(handles.asc_status(6) > 0)
    theta_KS_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(6)).Children(2).String);
end

% Criterion 7 has been chosen: Maximum deviation from NW OR
%'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
if(handles.asc_status(7) > 0)
    theta_NW_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(7)).Children(2).String);
end

%% get input for base vectors from GUI
base_aust = zeros(3,3);
k = 19; % position of 1st entry of 1st base-vec for austenite
for i = 1:3
    for j = 1:3
        base_aust(i,j) = str2num(handles.pan_base_vec.Children(k).String);
        k = k-1;
    end
end
austenite.my_base = base_aust;

base_mart = zeros(3,3);
k = 9; % position of 1st entry of 1st base-vec for martensite
for i = 1:3
    for j = 1:3
        base_mart(i,j) = str2num(handles.pan_base_vec.Children(k).String);
        k = k-1;
    end
end
martensite.my_base = base_mart;

                   
% get input for correspondence matrix from GUI
C_am = zeros(3,3);
k = 9; % counter for position in array handles.pan_corrmat.Children(k)
for i = 1:3
    for j = 1:3
        C_am(i,j) = str2num(handles.pan_corrmat.Children(k).String);
        k = k-1;
    end
end
martensite.C_am = C_am;

%%
% EHL: TODO: add functionality, to read input for slip planes and
% directions
plane_families =     [ [1 1 0] ;
                       [1 1 2] ];
direction_families = [ [1 1 1]; 
                       [1 1 0] ];


%% EHL: add input in GUI?
austenite.Bravais_type  = 'cubic';
martensite.Bravais_type = 'cubic';
% austenite.Centering = 
martensite.Centering = 'I';

austenite.Lp = a_aust*[1 1 1];  % 3.5975576 % {1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
                                              % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
martensite.Lp = a_mart *[1 1 1];  % 2.8807346  

                  
% further initializations
% define Bain-strain
eta1 = (a_mart/a_aust)*sqrt(2);
eta3 = a_mart / a_aust; % this is one form of three pcossible for the bain strain

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf 1.0 getuned
% werden. Die differenz ist also (n1-1).
B3 = [eta1 0    0   
       0  eta1  0
       0  0  eta3];
martensite.U = B3;

display('Volume change in percent is:');
detB3 = det(B3)

cpps_gamma = all_from_family_perms( [1 1 1] ); % close packed planes of gamma-lattice
KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
NW = all_from_family_perms( [1 2 1], false );