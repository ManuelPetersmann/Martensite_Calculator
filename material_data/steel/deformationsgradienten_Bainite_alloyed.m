function [ VariantenSort ] = deformationsgradienten_Bainite( ) 
% Formulas due to Bhadeshia Paper
% all Angels in Radians!!! Matlab Standard

clc
format long
% Investigate Paper: Davenport 1974

% 1. Aufstellen des Bain strains "U" (variable Gitterparameter)
% 2. Modifizieren von U durch Scherung "S" auf physikalisch sinnvollem Slip
% system
% 3. Berechnen von S*U 
% Festlegen der Drehachse auf die Habitusebene und überpfrüfen ob sich durch eine 
% Rotation "R" R*S*U auf die Angegebene Form eines invariant plane strains siehe
% Christian S 58 bringen lässt.

% lattice parameter of pure iron
a_fe = 2.8664 

% BCC Austenite constituents in weight percent
M_C = 0.1 %0.75
M_Si = 1.63
M_Mn = 1.95
M_Ni = 0.
M_Mo = 0.29
M_Cr = 0.1
M_V = 0.01
%
M_Al = 0
w_i = [M_C, M_Si, M_Mn, M_Ni, M_Mo, M_Cr, M_V];
size(w_i,2)

%%%% Molar Masses of constitutional elements in g/mol
MM_C = 12.011
MM_Si = 28.085
MM_Mn = 54.938
MM_Ni = 58.693
MM_Mo = 95.94
MM_Cr = 51.996
MM_V = 50.942
MM = [MM_C, MM_Si, MM_Mn, MM_Ni, MM_Mo, MM_Cr, MM_V];
x_i = zeros(1,size(w_i,2));
% convert weight percent in molar percent
M_total = sum(w_i./MM); % see wikipedia Mass-fraction --> Mole-fraction
for k = 1:size(w_i,2)
    x_i(k) = ( w_i(k) / MM(k) ) / M_total;
end

% calculate lattice parameters as functions of alloying elements
% 
a_bct = a_fe + ( (a_fe - 0.279*x_i(1) )^2 * (a_fe + 2.496*x_i(1) ) - a_fe^3 ) / ( 3*a_fe^2 )  ...
    - 0.03*x_i(2) + 0.06*x_i(3) + 0.07* x_i(4) + 0.31*x_i(5) + 0.05*x_i(6) + 0.096* x_i(7);
a_bct

c_bct = a_bct;

% 
a_fcc = 3.578 + 0.033*M_C + 0.00095*M_Mn - 0.0002*M_Ni + 0.0006*M_Cr + 0.0056*M_Al ...
    +0.0031*M_Mo + 0.0018*M_V;
a_fcc

% K-S orientation relationship: {111}a --> {110}m & <110>a --> <111>m
% observed for C > 0.5w
% N-W OR: {111}a --> {011}m & <112>a --> <011>m
% Unterschied dieser beiden OR sind 5,26° (Überführung durch rotation um
% richtige <011>m möglich. K-S is exactly midway between two N-W relations


% define Bain-streches
n1 = (a_bct/a_fcc)*sqrt(2));

%n2 = a_bct/(a_fcc/sqrt(2));

n2 = c_bct / a_fcc; % this is one form of three possible for the bain strain

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf eins getuned
% werden. Die differenz ist also (n1-1).

B = [n1 0  0   
     0  n1 0
     0  0  n2];
 
%  variant1 = B

det( B ) 

% Transformation matrix for the transformation from the principle axis
% system to the cubic system in order to determine all variants by applying
% all rotations to the deformation gradient that map the cube back to
% itshelf

% Correspondance Matrix
 C = [1/sqrt(2),  -1/sqrt(2), 0;
      1/sqrt(2),  1/sqrt(2),  0;
       0,       0,      1] ;

variant1 = C*B*C'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assamble all roationas that map the cubic lattice back to itself
Rot_Matrizen = cubic_pointgroup_rotations()

for i=1:size(Rot_Matrizen,3)
    Varianten(:,:,i) = (Rot_Matrizen(:,:,i)) * variant1 * Rot_Matrizen(:,:,i)';
end

VariantenSort = matrixArray_reduzieren(Varianten, 1e-6);

end


