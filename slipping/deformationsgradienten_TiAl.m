function [ VariantenSort, Rot_Matrizen ] = deformationsgradienten_TiAl( )
%all Angels in Radians!!! Matlab Standard
clc
format long

% lattice parameters of direct measurement at Ms of about 805�C
abeta  = 0.323; 
aalpha = 0.29; 
calpha = 0.4659;  
% former lattice parameters
% abeta  = 0.3277; 
% aalpha = 0.293; 
% calpha = 0.47;

n1 = (1/sqrt(2))*(calpha/abeta);
mittlerer_eigenwert = n1

n2 = (sqrt(3)/sqrt(2))*(aalpha/abeta);

n3 = (aalpha/abeta);

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf eins getuned
% werden. Die differenz ist also (n1-1).

B = [n1 0  0   
     0  n2 0
     0  0  n3];

trace( B ) 
 
%% Pitsch Schrader cartesian coordinate system relationship between crystallographic
% directions in bcc and hcp 5.26° 

% teta = 60-(180/pi)*acos(1/sqrt(3)) %Ergebnis acos in radians...
% teta = teta*(pi/180) % teta in radians 

% Rplus = [1      0        0
%         0   cos(teta)   -sin(teta)
%         0   sin(teta)   cos(teta)]
  
% Rminus = [1     0         0
%           0   cos(-teta)  -sin(-teta)
%           0   sin(-teta)   cos(-teta)]
%%

% Transformation matrix for the transformation from the principle axis
% system to the cubic system in order to determine all variants by applying
% all rotations to the deformation gradient that map the cube back to
% itshelf

%Bp1 = [1 0 1];
%Bp1 = Bp1 / norm( Bp1 );
%Bp2 = [0 1 0]; % Principal basis of deformation gradient
%Bp2 = Bp2 / norm( Bp2 );
%Bp = [Bp1; Bp2; cross(Bp1,Bp2)];
%Bc = [1 0 0; 0 1 0; 0 0 1]

%TT = Basiswechselmatrix( ( [1 0 1]/ norm( [1 0 1] )), ( [1 0 -1]/ norm( [1 0 -1] ) ), [1 0 0], [0 1 0] )

 T = [1/sqrt(2),  1/sqrt(2), 0;
      0,          0,         1;
      1/sqrt(2), - 1/sqrt(2), 0] 

TBTtrans = T*B*T'
trace( TBTtrans )

%T_B_Ttrans = T*B*T' %T*Rplus*B*T'
% T_Rminus_B_Ttrans = T*B*T' %T*Rminus*B*T'


%% in den nächsten Zeilen werden die 24 Rotationsmatrizen berechnet und in
% einem 3D array hintereinander abgespeichert

Rot_Matrizen = eye(3); % erster Eintrag ist die Identity, die 3. Dimension kommt im Nachhinein dazu

zaehler = 1;

    function[Rot] = rotieren(n, alpha)
alpha = alpha*(pi/180); %umrechnen von Grad aus Argument auf Radianten für cos und sin      
n1=n(1);
n2=n(2);
n3=n(3);
         
Rot = round( [(n1^2)*(1-cos(alpha))+cos(alpha)      n1*n2*(1-cos(alpha))-n3*sin(alpha)      n1*n3*(1-cos(alpha))+n2*sin(alpha)
        n2*n1*(1-cos(alpha))+n3*sin(alpha)   (n2^2)*(1-cos(alpha))+cos(alpha)         n2*n3*(1-cos(alpha))-n1*sin(alpha)
        n3*n1*(1-cos(alpha))-n2*sin(alpha)   n3*n2*(1-cos(alpha))+n1*sin(alpha)       (n3^2)*(1-cos(alpha))+cos(alpha)] ) ;
    end

display('die 9 Rotationen um die 3 Achsen')
e= [1 0 0 
    0 1 0
    0 0 1];

for i= 1:3
    alpha=90;
    for j = 1:3
        zaehler=zaehler+1;
        Rot_Matrizen(:,:,zaehler) = rotieren(e(:,i), alpha);        
        alpha=alpha+90;
    end
end

%display('Die 6 Rotationen um die Fl�chendiagonalen')
e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];
alpha=180;

n = 1/sqrt(2)*(e1+e2);
Rot_Matrizen(:,:,11)= rotieren(n, alpha);

n = 1/sqrt(2)*(e1-e2);
Rot_Matrizen(:,:,12)= rotieren(n, alpha);

n = 1/sqrt(2)*(e2+e3);
Rot_Matrizen(:,:,13)= rotieren(n, alpha);

n = 1/sqrt(2)*(e2-e3);
Rot_Matrizen(:,:,14)= rotieren(n, alpha);

n = 1/sqrt(2)*(e3+e1);
Rot_Matrizen(:,:,15)= rotieren(n, alpha);

n = 1/sqrt(2)*(e3-e1);
Rot_Matrizen(:,:,16)= rotieren(n, alpha);

%display('Die 8 Rotationen um die Raumdiagonalen')
alpha = 120;

n = 1/sqrt(3)*(e1+e2+e3);
Rot_Matrizen(:,:,17)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1+e2-e3);
Rot_Matrizen(:,:,18)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1-e2+e3);
Rot_Matrizen(:,:,19)= rotieren(n, alpha);

n = 1/sqrt(3)*(-e1+e2+e3);
Rot_Matrizen(:,:,20)= rotieren(n, alpha);

alpha=240;

n = 1/sqrt(3)*(e1+e2+e3);
Rot_Matrizen(:,:,21)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1+e2-e3);
Rot_Matrizen(:,:,22)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1-e2+e3);
Rot_Matrizen(:,:,23)= rotieren(n, alpha);

n = 1/sqrt(3)*(-e1+e2+e3);
Rot_Matrizen(:,:,24)= rotieren(n, alpha);

%display('Die 24 Rotationsmatrizen die sich aus der Punktsymmetrie ergeben')
Rot_Matrizen;
%%

%display('Die sich aus T_Rplus_B_Ttrans und T_Rminus_B_Ttrans ergebenen Varianten U mit U=sqrtm(F^t*F)')  

% U1 = sqrtm(T_Rplus_B_Ttrans' * T_Rplus_B_Ttrans);
% U25 = sqrtm(T_Rminus_B_Ttrans' * T_Rminus_B_Ttrans);

%Jetzt werden T_Rplus_B_Ttrans und T_Rminus_B_Ttrans mit allen 24
%Rotationsmatrizen rotiert. Auf diese Weise erh�lt man alle M�glichen
%Deformationsgradienten, in Summe 48, davon sind jeweils 4 gleich??? --�berpr�fen also 
%erhält man 12 verschiedene.
%Deformationsgradienten = eye(3)

for i=1:24
    %Deformationsgradienten(:,:,i) = (Rot_Matrizen(:,:,i))' * TBTtrans * Rot_Matrizen(:,:,i);   
    Varianten(:,:,i) = (Rot_Matrizen(:,:,i)) * TBTtrans * Rot_Matrizen(:,:,i)';
end
 
%display('Der vollständige Satz an Deformationsgradienten')
%Deformationsgradienten;
%display('Der vollständige Satz an Varianten')
%Varianten;

VariantenSort = matrixArray_reduzieren(Varianten, 1e-6);

%for i = 1:12
%VariantenvonDeformationsgradienten(:,:,i)= sqrtm((DeformationsgradientenSort(:,:,i))' * DeformationsgradientenSort(:,:,i) );
%end

%VariantenvonDeformationsgradientenSort = reduzieren(VariantenvonDeformationsgradienten);

%VariantenvonDeformationsgradientenSort;

%for i = 1:6
%det(VariantenSort(:,:,i))
%end

%for i = 1:6
%det(DeformationsgradientenSort(:,:,i))
%end
%{
F1 = DeformationsgradientenSort(:,:,2)

F2 = DeformationsgradientenSort(:,:,12)

zwischen1= F1'*F1

zwischen2 = F2'*F2

U1 = sqrtm(F1'*F1)

U2 = sqrtm(F2'*F2)
%}


end


