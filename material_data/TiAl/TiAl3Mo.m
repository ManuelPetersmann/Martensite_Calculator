% lattice parameters of direct measurement at Ms of about 805�C
abeta  = 0.323; 
aalpha = 0.29; 
calpha = 0.4659;  
% former lattice parameters
% abeta  = 0.3277; 
% aalpha = 0.293; 
% calpha = 0.47;

n1 = (1/sqrt(2))*(calpha/abeta);
mittlerer_eigenwert = n1;

n2 = (sqrt(3)/sqrt(2))*(aalpha/abeta);

n3 = (aalpha/abeta);

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf eins getuned
% werden. Die differenz ist also (n1-1).

B = [n1 0  0   % maybe this should be compared with structrans sometime...
     0  n2 0
     0  0  n3];

%eigs(B'*B)
% 1 - lambda_2 = 0.04 i.e. relative 4% 

 % Pitsch Schrader cartesian coordinate system relationship between crystallographic
% directions in bcc and hcp 5.26° 

teta = 60-(180/pi)*acos(1/sqrt(3)); %Ergebnis acos in radians...
teta = teta*(pi/180); % teta in radians 

Rplus = [1      0        0
        0   cos(teta)   -sin(teta)
        0   sin(teta)   cos(teta)];
  
Rminus = [1     0         0
          0   cos(-teta)  -sin(-teta)
          0   sin(-teta)   cos(-teta)];