function [ full, voigt ] = assembleElasticTensor( )

% From fourth order tensor to voight notation by
% 11-->1, 22-->2, 33-->3, 12-->4, 13-->5, 23-->6

% Voigt Matrix:
% [ 11, 12, 13, 14, 15, 16 ]
% [ 21, 22, 23, 24, 25, 26 ]
% [ 31, 32, 33, 34, 35, 36 ]
% [ 41, 42, 43, 44, 45, 46 ]
% [ 51, 52, 53, 54, 55, 56 ]
% [ 61, 62, 63, 64, 65, 66 ]
% this matrix is symmetric so there are at most 21 independent elastic
% constants

% Sehitoglu, Wang - Niti modulus dilemma
v11 = 141.e-9  % 2.54e-7
v12 = 74.e-9   % 1.04e-7
v22 = 179e-9   % 1.8e-7
v13 = 102.e-9  % 1.36e-7
v23 = 121e-9   % 1.51e-7
v33 = 184e-9   %2.48e-7
v14 = 0
v24 = 0
v34 = 0
v44 = 50.e-9   %9.1e-8
v15 = -43.e-9  %0
v25 = -24.e-9  %0
v35 = -1e-9    %0
v45 = 0.      %-3e-9
v55 = 30.e-9  %9.3e-8
v16 = 0.      %2.1e-8
v26 = 0.
v36 = 0.       %-6e-9
v46 = 1.5e-9   %0
v56 = 0
v66 = 67.e-9   %5e-9

% B19' constants thesis rostdam golesorktabar, coordinate system Quantum
% % espresso see thesis
% v11 = 2.54e-7
% v12 = 1.04e-7
% v22 = 1.8e-7
% v13 = 1.36e-7
% v23 = 1.51e-7
% v33 = 2.48e-7
% v14 = 0
% v24 = 0
% v34 = 0
% v44 = 9.1e-8
% v15 = 0
% v25 = 0
% v35 = 0
% v45 = -3e-9
% v55 = 9.3e-8
% v16 = 2.1e-8
% v26 = 0
% v36 = -6e-9
% v46 = 0
% v56 = 0
% v66 = 5e-9

% % Tiphaine Pelliset 4H-SiC, all values are GPa
% v11 = 501
% v12 = 111
% v22 = 501
% v13 = 52
% v23 = 52
% v33 = 553
% v14 = 0
% v24 = 0
% v34 = 0
% v44 = 163
% v15 = 0
% v25 = 0
% v35 = 0
% v45 = 0
% v55 = 163
% v16 = 0
% v26 = 0
% v36 = 0
% v46 = 0
% v56 = 0
% v66 = 195.5

% wagner windl
% v11 = 2.23e-7
% v12 = 1.29e-7
% v22 = 2.41e-7
% v13 = 9.9e-8
% v23 = 1.25e-7
% v33 = 2.00e-7
% v14 = 0
% v24 = 0
% v34 = 0
% v44 = 7.6e-8
% v15 = 2.7e-8
% v25 = -9.e-9
% v35 = 4.e-9
% v45 = 0
% v55 = 2.1e-8
% v16 = 0
% v26 = 0
% v36 = 0
% v46 = -4.e-9
% v56 = 0
% v66 = 7.7e-8

% v11 = 1.41e-7
% v12 = 7.4e-8
% v22 = 1.79e-7
% v13 = 1.02e-7
% v23 = 1.21e-7
% v33 = 1.84e-7
% v14 = 0
% v24 = 0
% v34 = 0
% v44 = 5.0e-8
% v15 = -4.3e-8
% v25 = -2.4e-8
% v35 = -1.e-9
% v45 = 0
% v55 = 3.0e-8
% v16 = 0
% v26 = 0
% v36 = 0
% v46 = 1.5e-9
% v56 = 0
% v66 = 6.7e-8

% v11 = input( 'Enter the value for c11 = C1111':    )
% v12 = input( 'Enter the value for c12 = C1122':    )
% v22 = input( 'Enter the value for c12 = C1122':    )
% v13 = input( 'Enter the value for c12 = C1122':    )
% v23 = input( 'Enter the value for c12 = C1122':    )
% v33 = input( 'Enter the value for c12 = C1122':    )
% v14 = input( 'Enter the value for c12 = C1122':    )
% v24 = input( 'Enter the value for c12 = C1122':    )
% v34 = input( 'Enter the value for c12 = C1122':    )
% v44 = input( 'Enter the value for c12 = C1122':    )
% v15 = input( 'Enter the value for c12 = C1122':    )
% v25 = input( 'Enter the value for c12 = C1122':    )
% v35 = input( 'Enter the value for c12 = C1122':    )
% v45 = input( 'Enter the value for c12 = C1122':    )
% v55 = input( 'Enter the value for c12 = C1122':    )
% v16 = input( 'Enter the value for c12 = C1122':    )
% v26 = input( 'Enter the value for c12 = C1122':    )
% v36 = input( 'Enter the value for c12 = C1122':    )
% v46 = input( 'Enter the value for c12 = C1122':    )
% v56 = input( 'Enter the value for c12 = C1122':    )
% v66 = input( 'Enter the value for c12 = C1122':    )


M = zeros(3,3,3,3);

% [ x1, x2, x3, x4, x5 ] = deal( 8 ); is equal to x1 = x2 = x3 = 8
M(1,1,1,1) = v11;
[ M(1,1,2,2), M(2,2,1,1)] = deal ( v12 );
M(2,2,2,2) = v22;
[ M(1,1,3,3), M(3,3,1,1) ] = deal ( v13 );
[ M(2,2,3,3) , M(3,3,2,2) ] = deal ( v23);
M(3,3,3,3) = v33;
[ M(1,1,1,2) , M(1,2,1,1) , M(1,1,2,1) , M(2,1,1,1) ] = deal( v14 ) ;
[ M(2,2,1,2) , M(1,2,2,2) , M(2,2,2,1) , M(2,1,2,2) ] = deal ( v24 ) ;
%second line
[ M(3,3,1,2) , M(1,2,3,3) , M(3,3,2,1) , M(2,1,3,3) ] = deal( v34 );
[ M(1,2,1,2) , M(2,1,1,2) , M(1,2,2,1) , M(2,1,2,1) ] = deal ( v44 );
[ M(1,1,1,3) , M(1,3,1,1) , M(1,1,3,1) , M(3,1,1,1) ] = deal( v15 );
[ M(2,2,1,3) , M(1,3,2,2) , M(2,2,3,1) , M(3,1,2,2) ] = deal( v25 );
[ M(3,3,1,3) , M(1,3,3,3) , M(3,3,3,1) , M(3,1,3,3) ] = deal( v35 );
[ M(1,2,1,3) , M(1,3,1,2) , M(2,1,1,3) , M(1,2,3,1) , M(3,1,1,2) , M(1,3,2,1) , M(2,1,3,1) , M(3,1,2,1) ] = deal( v45 );
[ M(1,3,1,3) , M(3,1,1,3) , M(1,3,3,1) , M(3,1,3,1) ] = deal( v55);
[ M(1,1,2,3) , M(2,3,1,1) , M(1,1,3,2) , M(3,2,1,1) ] = deal( v16 );
%third line
[ M(2,2,2,3) , M(2,3,2,2) , M(2,2,3,2) , M(3,2,2,2) ] = deal( v26 );
[ M(3,3,2,3) , M(2,3,3,3) , M(3,3,3,2) , M(3,2,3,3) ] = deal( v36 );
[ M(1,2,2,3) , M(2,3,1,2) , M(2,1,2,3) , M(1,2,3,2) , M(3,2,1,2) , M(2,3,2,1) , M(2,1,3,2) , M(3,2,2,1) ] = deal( v46 );
[ M(1,3,2,3) , M(2,3,1,3) , M(3,1,2,3) , M(1,3,3,2) , M(3,2,1,3) , M(2,3,3,1) , M(3,1,3,2) , M(3,2,3,1) ] = deal( v56 );
[ M(2,3,2,3) , M(3,2,2,3) , M(2,3,3,2) , M(3,2,3,2) ] = deal( v66 );

full = M;

voigt = [ M(1,1,1,1) , M(1,1,2,2) , M(1,1,3,3) , M(1,1,1,2) , M(1,1,1,3) , M(1,1,2,3) ;
                M(2,2,1,1) , M(2,2,2,2) , M(2,2,3,3) , M(2,2,1,2) , M(2,2,1,3) , M(2,2,2,3) ;
                M(3,3,1,1) , M(3,3,2,2) , M(3,3,3,3) , M(3,3,1,2) , M(3,3,1,3) , M(3,3,2,3) ;
                M(1,2,1,1) , M(1,2,2,2) , M(1,2,3,3) , M(1,2,1,2) , M(1,2,1,3) , M(1,2,2,3) ;
                M(1,3,1,1) , M(1,3,2,2) , M(1,3,3,3) , M(1,3,1,2) , M(1,3,1,3) , M(1,3,2,3) ;
                M(2,3,1,1) , M(2,3,2,2) , M(2,3,3,3) , M(2,3,1,2) , M(2,3,1,3) , M(2,3,2,3) ];

% C = zeros(3,3,3,3)
% 
% voigt_entries = [v11, v12, v22, v13, v23, v33, v14, v24, v34, v44, v15, v25, v35, v45, v55, v16, v26, v36, v46, v56, v66]

% voight_indices is a tuple with the entry and a list of indixes at each position
% voigt_indices = [ 1,1,1,1;
%                   1,1,2,2;
%                   2,2,2,2;
%                   1,1,3,3;
%                   2,2,3,3;
%                   3,3,3,3;
%                   1,1,1,2;
%                   2,2,1,2;
%                   3,3,1,2;
%                   1,2,1,2;
%                   1,1,1,3;
%                   2,2,1,3;
%                   3,3,1,3;
%                   1,2,1,3;
%                   1,3,1,3;
%                   1,1,2,3;
%                   2,2,2,3;
%                   3,3,2,3;
%                   1,2,2,3;
%                   1,3,2,3;
%                   2,3,2,3];

% Operationen:
% Due to the symmetry of stress and strain the number of independent
% entries reduces to 36: Cijkl = Cjikl = Cijlk = Cjilk
% Moreover, the existence of a unique strain energy potential requires that
% Cijkl = Cklij
% Gegeben ist einer der 21 verschiedenen Eintraege der Compliance
% 1 Vertauschen des ersten und zweiten Paares (Cijkl = Cjikl % Spannungstensorsymmetrie) --> 54 entries)
% 2 Vertauschen des ersten und zweiten index (Cijkl = Cijlk Dehnungstensorsymmetrie --> 36 entries)
% 3 Achtung es gehört auch die Kombinierte operation dazu also 1 + 2 Cijkl = Cjilk
% 4 The strain energy function should not change when we interchange ij and kl in the quadratic form, we must have Cijkl = Cklij
              
% for count = 1 : 21
% 		% give indizes variables
% 		i = voigt_indices(count, 1);
%         j = voigt_indices(count, 2);
%         k = voigt_indices(count, 3);
%         l = voigt_indices(count, 4);
% 		% set given entry on right place of tensor
% 		entry = voigt_entries(count);
% 		% carry out 6 other symmetry operations
% 		% 1 swap index pairs:
% 		if C(k,l,i,j) == 0
% 			C(k,l,i,j) = entry;
%         end
% 		% 2 swap first two first indizes
% 		if C(j,i,k,l) == 0
%             C(j,i,k,l) = entry;
%         end
%         % 3 swap last two indizes
% 		if C(i,j,l,k) == 0
% 			C(i,j,l,k) = entry;
%         end
% 		% 4 = 1+2
% 		if C(l,k,i,j) == 0
%             C(l,k,i,j) = entry;
%         end
% 		% 5 = 1+3
% 		if C(k,l,j,i) == 0
% 			C(k,l,j,i) = entry;
%         end
% 		% 6 = 2+3
% 		if C(j,i,l,k) == 0
% 			C(j,i,l,k) = entry;
%         end
% 		% 7 = 1+2+3
% 		if C(l,k,j,i) == 0
% 			C(l,k,j,i) = entry;
%         end

end

