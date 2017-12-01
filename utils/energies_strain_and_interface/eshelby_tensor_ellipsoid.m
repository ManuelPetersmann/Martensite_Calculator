function S = eshelby_tensor_ellipsoid( nu, a,b,c )
% call: eshelby_tensor_ellipsoid(G,nu,a,b,c)
% This function takes the Poissons ratio....nu and the half axis lengths
% of an ellipsoid: a > b > c and calculates the Eshelby tensor

a2 = a^2;
b2 = b^2;
c2 = c^2;
l2 = [a2,b2,c2];

theta = 1./sin( sqrt(1.- c2/a2) );

k = sqrt(a2 - b2) / sqrt(a2 - c2);

% solve incomplete elliptic integral of the second kind E(theta,k) 
% https://de.mathworks.com/help/symbolic/mupad_ref/elliptice.html
E = ellipticE(theta,k);

% solve incomplete elliptic integral of the first kind F(theta,k) 
% https://de.mathworks.com/help/symbolic/mupad_ref/ellipticf.html
F = ellipticF(theta,k);

pre = 4*pi*a*b*c;

Q = 3./(8*pi*(1.-nu));

R = (1.-2*nu) / (8*pi*(1.-nu));

%% ellipsoid
% Ia = pre * (F-E) / ( (a2-b2)*sqrt(a2-c2) );
% 
% Ic = pre / ( (b2-c2)*sqrt(a2-c2) ) * ( b*sqrt(a2-c2) / (a*c) - E);
% 
% Ib = 4*pi - Ia - Ic;

% I_i = [Ia, Ib, Ic];
% 
% % I_ij = [ Iab, Iac, Iba, Ibc, Ica, Icb ]
% % indizes:   1    2    3    4    5    6
% index = 0;
% for i = 1:3
%     for j = 1:3
%         if abs(i-j) > 1.e-3 % if i-not-j
%             index = index + 1;
%             I_ij(index) = (I_i(j) - I_i(i) ) / (3*( l2(i) - l2(j) ) );
%         end
%     end
% end
% 
% I_kk = [I_aa, I_bb, I_cc]
% pre2 = 3.*pi / (3.*a2);
% for k = 1:3
%    I_kk(k) = pre2 - I_ij( 2*k-1 ) - I_ij( 2*k );
% end
% 
% S = zeros(3,3,3,3);

%% elliptic cylinder after original paper eshelby 1957
Ia = 4*pi*b / (a+b);
Ib = 4*pi*a / (a+b);
Ic = 0;
I_ab = 4*pi/(3*(a+b)); % not sure in paper - could be 4*pi/3*(a+b)...
%
I_aa = 4*pi/(3*a^2) - I_ab; % not sure - could be 4*pi/(3*a^2 - I_ab);
I_bb = 4*pi/(3*b^2) - I_ab; % deto
I_ac = 0.;
I_bc = 0.;
I_cc = 0.;
I_kk = [I_aa, I_bb, I_cc];
I_ij = [ I_ab, I_ac, 0, I_bc, 0, 0 ];

%%

% [ x1, x2, x3, x4, x5 ] = deal( 8 ); is equal to x1 = x2 = x3 = 8
%
%
S(1,1,1,1) = Q*a2*I_kk(1) + R*Ia;
S(2,2,2,2) = Q*b2*I_kk(2) + R*Ib;
S(3,3,3,3) = Q*c2*I_kk(3) + R*Ic;
% relations of the form S1122 = S2211 are not valid - check!
S(1,1,2,2) = Q* b2* I_ij(1) - R*Ia;
S(1,1,3,3) = Q* c2* I_ij(2) - R*Ia;
S(2,2,1,1) = Q* a2* I_ij(3) - R*Ib;
S(2,2,3,3) = Q* c2* I_ij(4) - R*Ib;
S(3,3,1,1) = Q* a2* I_ij(5) - R*Ic;
S(3,3,2,2) = Q* b2* I_ij(6) - R*Ic;
%
% The Eshelby tensor satisfies minor symmetries
% S_ijkl = S_jikl = S_ijlk   
% but in general no major symmetries
% S_ijkl not S_klij - see e.g. Bower solidmechanics
%
[ S(1,2,1,2) , S(2,1,1,2) ] = deal( 0.5*Q*(a2+b2)*I_ij(1) + 0.5*R*(Ia + Ib) );
% S(1,2,2,1), S(2,1,2,1) - would have major symmetry  

[ S(2,3,2,3) , S(3,2,2,3) ] = deal( 0.5*Q*(b2+c2)*I_ij(4) + 0.5*R*(Ib + Ic) );
% S(2,3,3,2) , S(3,2,3,2) - would have major symmetry

[ S(1,3,1,3) , S(3,1,1,3) ] = deal( 0.5*Q*(c2+a2)*I_ij(5) + 0.5*R*(Ic + Ia) );
% S(1,3,3,1) , S(3,1,3,1) - would have major symmetry
%
% Coefficients coupling an extension and a shear (S1112, S1123, S2311 ...)
% or one shear to another (S1223...) are zero. However same shears S1212 not 0! 
%
%
%
%
%[ S(1,1,1,2) , S(1,2,1,1) , S(1,1,2,1) , S(2,1,1,1) ] = deal( 0. ) ;
%[ S(2,2,1,2) , S(1,2,2,2) , S(2,2,2,1) , S(2,1,2,2) ] = deal ( 0. ) ;
%second line
%[ S(3,3,1,2) , S(1,2,3,3) , S(3,3,2,1) , S(2,1,3,3) ] = deal( 0. );
%[ S(1,2,1,2) , S(2,1,1,2) , S(1,2,2,1) , S(2,1,2,1) ] = deal( 0. );
%[ S(1,1,1,3) , S(1,3,1,1) , S(1,1,3,1) , S(3,1,1,1) ] = deal( 0. );
%[ S(2,2,1,3) , S(1,3,2,2) , S(2,2,3,1) , S(3,1,2,2) ] = deal( 0. );
%[ S(3,3,1,3) , S(1,3,3,3) , S(3,3,3,1) , S(3,1,3,3) ] = deal( 0. );
%[ S(1,2,1,3) , S(1,3,1,2) , S(2,1,1,3) , S(1,2,3,1) , S(3,1,1,2) , S(1,3,2,1) , S(2,1,3,1) , S(3,1,2,1) ] = deal( 0. );
%[ S(1,3,1,3) , S(3,1,1,3) , S(1,3,3,1) , S(3,1,3,1) ] = deal( 0.);
%[ S(1,1,2,3) , S(2,3,1,1) , S(1,1,3,2) , S(3,2,1,1) ] = deal( 0.);
%third line
%[ S(2,2,2,3) , S(2,3,2,2) , S(2,2,3,2) , S(3,2,2,2) ] = deal( 0. );
%[ S(3,3,2,3) , S(2,3,3,3) , S(3,3,3,2) , S(3,2,3,3) ] = deal( 0. );
%[ S(1,2,2,3) , S(2,3,1,2) , S(2,1,2,3) , S(1,2,3,2) , S(3,2,1,2) , S(2,3,2,1) , S(2,1,3,2) , S(3,2,2,1) ] = deal( 0. );
%[ S(1,3,2,3) , S(2,3,1,3) , S(3,1,2,3) , S(1,3,3,2) , S(3,2,1,3) , S(2,3,3,1) , S(3,1,3,2) , S(3,2,3,1) ] = deal( 0. );
%[ S(2,3,2,3) , S(3,2,2,3) , S(2,3,3,2) , S(3,2,3,2) ] = deal( 0. );

%full = S;

voigt = [ S(1,1,1,1) , S(1,1,2,2) , S(1,1,3,3) , S(1,1,1,2) , S(1,1,1,3) , S(1,1,2,3) ;
          S(2,2,1,1) , S(2,2,2,2) , S(2,2,3,3) , S(2,2,1,2) , S(2,2,1,3) , S(2,2,2,3) ;
          S(3,3,1,1) , S(3,3,2,2) , S(3,3,3,3) , S(3,3,1,2) , S(3,3,1,3) , S(3,3,2,3) ;
          S(1,2,1,1) , S(1,2,2,2) , S(1,2,3,3) , S(1,2,1,2) , S(1,2,1,3) , S(1,2,2,3) ;
          S(1,3,1,1) , S(1,3,2,2) , S(1,3,3,3) , S(1,3,1,2) , S(1,3,1,3) , S(1,3,2,3) ;
          S(2,3,1,1) , S(2,3,2,2) , S(2,3,3,3) , S(2,3,1,2) , S(2,3,1,3) , S(2,3,2,3) ];


end

