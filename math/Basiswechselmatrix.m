function [ T ] = Basiswechselmatrix( o1, o2, n1, n2 )
% given two cartesian coordinate systems defined by the basis set vectors
% o_i (old system) and n_i (new system) computes the transformation matrix
% for similarity transformations (only!!!) that transforms the old into the
% new one.
% generally the new coordinate system must be presented in the old
% coordinate system via the direction cosines form the transformation
% matrix

% check if its a mutually perpendicular cartesian coordinate system
% since also random numbers are used its checked against an error
% threshhold < 1e-6
if dot(n1,n2) > 1e-3 || dot(o1,o2) > 1e-3
     error('manuel_checkCartesian', 'the specified coordinate axis are not mutally perpendicular')
end

% define the third vector that defines the cartesian coordinate system
% and write the old and new basis into matrices
n3 = cross(n1,n2);
o3 = cross(o1,o2);
N = cat(1,n1,n2,n3);
O = cat(1,o1,o2,o3);

T = zeros(3,3);

for i = 1:3
    for j = 1:3
        T(i,j) = dot( N(i,:) , O(j,:) ) / ( norm( N(i,:) ) * norm( O(j,:) ) ); % definition of direction cosines!
    end
end

% Mathematicians call the transposed of T the Übergangsmatrix
% P = T';
% and then define the transformation generally as
% Transformedmatrix = P^-1 * initialMatrix * P
% Note that for my transformation matrix R'*R = R*R' = I holds!

% check if the transformation matrix is properly defined
% the determinand has to be 1
if det(T)-1 > 1e-5
    error('myApp:basiswechel', 'The determinant of the Übergangsmatrix is unequal to 1')
end
% Its transpose and its inverse must be the same
if inv(T) - T.' < 1e-5
    ;
else
    error('myApp:basiswechel', 'The inverse of the Übergangsmatrix and its transposed are not the same')
end
% if the standartbasis is B' then T = A ( here A'). see wikipedia
% Basiswechsel %A' if the standartbasis is B then T = inv(B')

% Diese Matrix ist quadratisch und invertierbar.
% Ihre Inverse beschreibt den Basiswechsel von B' zurück nach B.

% Probe - es sollte T = T2 sein! siehe Wikipedia Basiswechsel
B(:,1) = o1;
B(:,2) = o2;
B(:,3) = o3;
B_strich(:,1) = n1;
B_strich(:,2) = n2;
B_strich(:,3) = n3;

%T2 = inv(B_strich) * B % so lässt es sich ja noch einfacher berechnen :-)
% Insbesondere gilt: Ist B die Standardbasis, so gilt T = B'^-1
% Standartbasis ist z.B. [1 0 0], [0 1 0], [0 0 1]
end

