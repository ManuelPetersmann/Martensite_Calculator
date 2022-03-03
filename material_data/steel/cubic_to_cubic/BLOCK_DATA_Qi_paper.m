function [ VariantenSort, strainSort  ] = BLOCK_DATA_Qi_paper( )  %[ blockpairs, n1, n2]
%all Angels in Radians!!! Matlab Standard

clc
format short

% Qi, Khachaturyan Paper 2014.
% Relative Änderung Volumsdehnungen:
 a_aust = 3.56 % see Qi et al. 2014
 a_mart = 2.869 % see Qi et al. 2014
n1 = (a_mart/a_aust)*sqrt(2);
n3 = a_mart / a_aust; % this is one form of three possible for the bain strain

% temprary: rounding of values for direct comparison with results from paper
% PET: 7.11.17 - rounding of values leads to 0.49% change which is even more!
% n1 = round ( ((a_mart/a_aust)*sqrt(2))*1000 )/1000;
% n3 = round ( (a_mart / a_aust)*1000 )/1000;

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf eins getuned
% werden. Die differenz ist also (n1-1).
B3_Qi = [n1 0    0   
       0  n1  0
       0  0  n3];
GL_B3_qi = 0.5*(B3_Qi'*B3_Qi) - eye(3)
% 
 display('volume change Bain - Qi')
 det(B3_Qi)
 
%% 
% A1 = [1.0928   0.1034  0.1022;   
%       0.0928   1.1034  0.1022;  
%      -0.1369  -0.1581  0.8464];
%  
%  eigs(A1)
%  
% display('volume change - shape transformation variant1 - Qi')
% det(A1)

%%  
block = [1.0981   0.0981  0.1022;   
         0.0981   1.0981  0.1022;  
        -0.1475  -0.1475  0.8464];
%     
%[W1,D1] = eigs(block)
%det(W1)

%[~,R] = polardecomposition( block );
%[vec4(4), vec4(1:3)]  = rotmat_to_axis_angle(R); %vrrotmat2vec( R );
%vec4
block = block * 1.00135 % tune determinant to one
E = 0.5*(block'*block) - eye(3)
frob(E)
%[W2,D2] = eigs(E)

%det(W1)
%det(W2)

%RFE = W1*inverse(W2);
%[vec4RFE(4), vec4RFE(1:3)]  = rotmat_to_axis_angle(RFE); %vrrotmat2vec( R );
%vecRFE4

display('volume change - shape transformation variant block - Qi')
det(block)
det(block * 1.00135)

( det(B3_Qi) - det(block) )*100

block_modified_5thvalue = [1.09819   0.09819  0.10229;   
                           0.09819   1.09819  0.10229;  
                          -0.14759  -0.14759  0.84649];
 
dbmod = det(block_modified_5thvalue)
error_dV = ( det(B3_Qi) - dbmod ) %*100 

nu = 0.33
G = 210.e3;
K = compression_modulus( nu, G, 'nuG')
SENER = K*error_dV
%A2 = (block - eye(3) )*0.5;
%A1 = A2 +eye(3);
%det(A1)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble all rotations that map the cubic lattice back to itself
% Rot_Matrizen = cubic_pointgroup_rotations();
% 
% for i=1:size(Rot_Matrizen,3)
%     block_variants(:,:,i) = (Rot_Matrizen(:,:,i)) * block * Rot_Matrizen(:,:,i)';
% end
% 
% VariantenSort = matrixArray_reduzieren(block_variants, 1e-4);
% 
% % it turned out that i get all 24 Block variants, and I found that 
% % I started with shape deformation A_1^1 Eq. 38 as in Paper and my 12 is his second block variant
% % J.W. Morris Paper, dislocated microstructure...
% operator = VariantenSort(:,:,12)*inv(VariantenSort(:,:,1));
% for i=1:size(Rot_Matrizen,3)
%     all_operators(:,:,i) = (Rot_Matrizen(:,:,i)) * operator * Rot_Matrizen(:,:,i)';
% end
% 
% blockpairs = zeros(12,2);
% count = 1;
% for i = 1:24
%     for j = i:24
%         for k = 1:size(all_operators,3)
%             if VariantenSort(:,:,i) - all_operators(:,:,k) * VariantenSort(:,:,j) < 1.e-6
%                 blockpairs( count, : ) = [i,j];
%                 count = count +1;
%             end
%         end
%     end
% end
% 
% blockpairs = unique(blockpairs,'rows');
% 
% % ave = zeros(3,3,12);
% % epsilon = zeros(12,1);
% % a1 = zeros(3,12);
% % a2 = zeros(3,12);
% % n1 = zeros(3,12);
% % n2 = zeros(3,12);
% % Q1 = zeros(3,3,12);
% % Q2 = zeros(3,3,12);
% % 
% % for i=1:12
% %     ave(:,:,i) = 0.5 *(  VariantenSort(:,:,blockpairs(i,1)) + VariantenSort(:,:,blockpairs(i,2)) )
% %     [epsilon(i), a1(i), a2(i), n1(i), n2(i), Q1(:,:,i), Q2(:,:,i) ] = rank_one_kachaturyan( ave(:,:,i) );
% % end
% 
% fileID = fopen('block_data_Qi-Paper.txt', 'a');
% format_block_variant = '%7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \t %7.4f \n'; %\t %7.4f \t %7.4f \t %7.4f \t \n';
% for i=1:12
%     % calculate strains (infinite deformation)
%     strainSort(:,:,blockpairs(i,1)) = 0.5*( VariantenSort(:,:,blockpairs(i,1))' + VariantenSort(:,:,blockpairs(i,1)) ) - eye(3) ;
%     strainSort(:,:,blockpairs(i,2)) = 0.5*( VariantenSort(:,:,blockpairs(i,2))' + VariantenSort(:,:,blockpairs(i,2)) ) - eye(3) ;
%     fprintf(fileID,format_block_variant, get_order( strainSort(:,:,blockpairs(i,1)) ) );
%     fprintf(fileID,format_block_variant, get_order( strainSort(:,:,blockpairs(i,2)) ) );
% end
% fclose('all');
% 
%     function six_entries = get_order( M )
%         six_entries = [M(1,1) M(2,2) M(3,3) M(1,2) M(2,3) M(1,3)]; % Order definition as in my Zebfront Code % M(2,1) M(3,2) M(3,1)];
%     end
% 
% end

