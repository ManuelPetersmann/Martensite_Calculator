function S = eshelby_tensor_cylinder( nu, a,b )
% call: eshelby_tensor_cylinder(G,nu,a,b,c)

%amnu = 1.-nu;
ab = a + b;
znu1 = 1. / (2*nu - 1.);
am2nu = 1.-2*nu;
a2 = a^2;
b2 = b^2;

% formulas from mechanical behavior of materials - zaoui, pineau
S(1,1,1,1) = znu1 * (  (b2+2*a*b) / (ab^2)  +  b*am2nu / ab );
S(2,2,2,2) = znu1 * (  (a2+2*a*b) / (ab^2)  +  a*am2nu / ab );
S(3,3,3,3) = 0.;
S(1,1,2,2) = znu1 * (  b2 / ab^2   +   b*am2nu / ab );
S(2,2,1,1) = znu1 * (  a2 / ab^2   +   a*am2nu / ab );

S(1,1,3,3) = znu1 * 2*b*nu / ab;
S(3,3,1,1) = 0.;
S(2,2,3,3) = znu1 * 2*a*nu / ab;
S(3,3,2,2) = 0.;

% The Eshelby tensor satisfies minor symmetries
% S_ijkl = S_jikl = S_ijlk = S_jilk   
% but in general no major symmetries
% S_ijkl not S_klij - see e.g. Bower solidmechanics
[ S(1,2,1,2), S(2,1,1,2)] = deal ( znu1 * ( (a2+b2) / (2*ab^2) + am2nu/2.));
% S(1,2,2,1), S(2,1,2,1) - would have major symmetry  

[ S(1,3,1,3) , S(3,1,1,3) ] = deal( b / (2*ab) );
% S(1,3,3,1) , S(3,1,3,1) - would have major symmetry

[ S(2,3,2,3) , S(3,2,2,3) ] = deal( a / (2*ab) );
% S(2,3,3,2) , S(3,2,3,2) - would have major symmetry

%% other terms are zero per default
% Coefficients coupling an extension and a shear (S1112, S1123, S2311 ...)
% or one shear to another (S1223...) are zero. However same shears S1212 
% not 0! 

% [ S(1,1,1,2) , S(1,2,1,1) , S(1,1,2,1) , S(2,1,1,1) ] = deal( 0. ) ;
% [ S(2,2,1,2) , S(1,2,2,2) , S(2,2,2,1) , S(2,1,2,2) ] = deal ( 0. ) ;
% %second line
% [ S(3,3,1,2) , S(1,2,3,3) , S(3,3,2,1) , S(2,1,3,3) ] = deal( 0. );
% [ S(1,2,1,2) , S(2,1,1,2) , S(1,2,2,1) , S(2,1,2,1) ] = deal( 0. );
% [ S(1,1,1,3) , S(1,3,1,1) , S(1,1,3,1) , S(3,1,1,1) ] = deal( 0. );
% [ S(2,2,1,3) , S(1,3,2,2) , S(2,2,3,1) , S(3,1,2,2) ] = deal( 0. );
% [ S(3,3,1,3) , S(1,3,3,3) , S(3,3,3,1) , S(3,1,3,3) ] = deal( 0. );
% [ S(1,2,1,3) , S(1,3,1,2) , S(2,1,1,3) , S(1,2,3,1) , S(3,1,1,2) , S(1,3,2,1) , S(2,1,3,1) , S(3,1,2,1) ] = deal( 0. );
% [ S(1,3,1,3) , S(3,1,1,3) , S(1,3,3,1) , S(3,1,3,1) ] = deal( 0.);
% [ S(1,1,2,3) , S(2,3,1,1) , S(1,1,3,2) , S(3,2,1,1) ] = deal( 0.);
% %third line
% [ S(2,2,2,3) , S(2,3,2,2) , S(2,2,3,2) , S(3,2,2,2) ] = deal( 0. );
% [ S(3,3,2,3) , S(2,3,3,3) , S(3,3,3,2) , S(3,2,3,3) ] = deal( 0. );
% [ S(1,2,2,3) , S(2,3,1,2) , S(2,1,2,3) , S(1,2,3,2) , S(3,2,1,2) , S(2,3,2,1) , S(2,1,3,2) , S(3,2,2,1) ] = deal( 0. );
% [ S(1,3,2,3) , S(2,3,1,3) , S(3,1,2,3) , S(1,3,3,2) , S(3,2,1,3) , S(2,3,3,1) , S(3,1,3,2) , S(3,2,3,1) ] = deal( 0. );
% [ S(2,3,2,3) , S(3,2,2,3) , S(2,3,3,2) , S(3,2,3,2) ] = deal( 0. );

%full = S;

voigt = [ S(1,1,1,1) , S(1,1,2,2) , S(1,1,3,3) , S(1,1,1,2) , S(1,1,1,3) , S(1,1,2,3) ;
          S(2,2,1,1) , S(2,2,2,2) , S(2,2,3,3) , S(2,2,1,2) , S(2,2,1,3) , S(2,2,2,3) ;
          S(3,3,1,1) , S(3,3,2,2) , S(3,3,3,3) , S(3,3,1,2) , S(3,3,1,3) , S(3,3,2,3) ;
          S(1,2,1,1) , S(1,2,2,2) , S(1,2,3,3) , S(1,2,1,2) , S(1,2,1,3) , S(1,2,2,3) ;
          S(1,3,1,1) , S(1,3,2,2) , S(1,3,3,3) , S(1,3,1,2) , S(1,3,1,3) , S(1,3,2,3) ;
          S(2,3,1,1) , S(2,3,2,2) , S(2,3,3,3) , S(2,3,1,2) , S(2,3,1,3) , S(2,3,2,3) ];


end

