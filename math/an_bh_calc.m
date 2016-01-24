function [a,n] = an_bh_calc(y1, y3, e1, e3, X, Ui)
% given a matrix with eigenvalues y1 < 1, y2 = 1, y3 > 1
% and corresponding eigenvectors e1, e2, e3, 
% calculated Algorithm to calculate Twin plane shear and normal a and n
% as well as korresponding values b and h for the habit plane
% X = Tschi, r = ro, n = n mit Dach
%If Ct ~= (not equal) to I, shear a and normal n are given by:
%
%Aufruf: an_bh_calc(y1, y3, e1, e3, X, Ui)
% y1, y3 = kleinster und größter Eigenwert der Ct bzw Ch Matrix
% e1, e3 = die zu den Eigenwerten korrespondierenden Eigenvektoren
% X = Vorzeichen, dass aus rechnung mitgenommen wird
% Ui = die erste Variantenmatrix ist nur bei der Berechnung auf Twin Ebene
% hinzuzufügen
a1 = ((y3*(1-y1))/(y3-y1))^(1/2)*e1; % die nächsten 5: Zeilenvektoren
%
a3 = X*(y1*(y3-1)/(y3-y1))^(1/2)*e3;
% 
nhilf1 = (-(1-y1)^(1/2))*e1;
%
nhilf3 = X*((y3-1)^(1/2))*e3;
%
nhilfsvektor = nhilf1 + nhilf3;
%
if nargin == 6
    n_unnormiert = ((y3^(1/2)-y1^(1/2))/((y3-y1)^(1/2)))*Ui*nhilfsvektor;
else
    n_unnormiert = ((y3^(1/2)-y1^(1/2))/((y3-y1)^(1/2)))*nhilfsvektor;
end
%
roh = norm(n_unnormiert); % norm(n) = norm(n,2) Euklidische Norm
%
n = n_unnormiert/roh;
%
% X...+/- 1, r unequal 0, konstante such that |n| =1 and e1,e3 normalized
% eigenvectors of C1 corresponding to y1, y3.
a_unnormiert = (a1 + a3);
%
a = roh*a_unnormiert;
%
end

