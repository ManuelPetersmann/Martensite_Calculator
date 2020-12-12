function [ min_g_U_j_out, WW_out ] = trim_eigenvalue_to_unity( B, m, l, variantnr, min_g_U_j_in, WW_in )

format long

fileID = fopen('Habitplane_evaluation.txt', 'a');

if variantnr == 1;
fprintf(fileID,'Subsequently the 10 calculated slip systems with the lowest shear necessary \n');
fprintf(fileID, 'to obtain an eigenvalue of the matrix C_ij (see Ball and James) \n');
fprintf(fileID, 'necessary for an invariant plane, for each martensite variant are shown \n');
fprintf(fileID, 'm...normal to the shear plane, l...shearing direction, g...simple shear coefficient \n \n');
end

% U_j = zeros(3) % U_j must be defined in the workspace outside the
% function call to keep it in the memory between function calls!
shear_amount = [0, 0];
eigenvalues = [0 0 0];
nr = 1; % counting nr to save all solutions 2 for each slip system {m,l}
%eigenvalues = sort( eigs( B ) );

 for i = 1:size(m,1) % since there are two solutions for g for each slip system respectively
     m(i,:);
     l(i,:);
     % m-Inequality: Wenn ich m_inequality zum aussortieren mitnehme
     % bekomme ich erstens im Gegensatz zum Buch heraus, dass das Kriterium
     % nicht erfÃ¼llt ist und zweitens werden dann meine g=0 und werden
     % somit als kleinste g's identifiziert. also lasse ich es weg
     %      if m_inequality( eigenvalues, m(i,:) ) > 0
     %          display('m-criterion not fullfilled')
     %          continue
     %      end

     %% calculate the representation of the Bain strain Tensor in the shear
     % coordinate system. Therefore get the transformation matrix that
     % transforms the specified U_i into the basis of the shear system:
     R = Basiswechselmatrix( [1 0 0], [0 1 0], ( l(i,:)/ norm( l(i,:) ) ), ( m(i,:)/ norm( m(i,:) ) ) );
     U_rot_m_l = R * B * R';
     trace( U_rot_m_l );
     
     % the eigenvalues must be invariant on a coordinate transformation (check here)
     if eigs( U_rot_m_l ) - eigs( B ) > 1e-8
         error( 'myApp:basiswechel', 'the eigenvalues of the linear transformation matrix are not the same after the similarity transformation')
     end
     % Another check is to verify that the trace of the transformation is
     % invarinat which must be fullfilled
     if trace( U_rot_m_l) - trace( B ) > 1e-8
         error( 'myApp:basiswechel', 'the traces of the linear transformation matrix before and after the transformation are not equal')
     end
     
     %% define a symbolic matrix representing the simple shear in the coordinate
     % system m, l, k of the shear system. If m is the normal to the plane
     % (z or 3 axis) and l is the shear direction (x or 1) axis, the shear
     % matrix writes in the general form:
     syms g;
     G = sym( [1 0 g; 0 1 0; 0 0 1] ) ;
     
     % Multiply the additional simple shear onto the Bain strain to trim one of
     % the eigenvalues to unity
     F = sym( U_rot_m_l * G );      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % Note that for the singular-value function the following constructed
     % symmetric matrix is used. Note that this is the same as F*F = F' * F
     % = F * F' See also the condition of the matrix C_ij at Ball and James
     % C_ij = U_i^-1 * U_j^2 * U_i^-1. Daher muss diese Gleichung verwendet
     % werden. Dies ist genau die Singulärwertgleichung und die
     % Singulärwerte haben im Gegensatz zu den Eigenwerten die Geometrische
     % Deutung die Halbachsen der verzerrten Elipse zu sein C = F^2; % Note
     % that .' is the transpose and only ' is the conjugate transpose!!!
     
     % Die erste Lösung wenn man C = F^-2 statt F nimmt stimmt mit der von
     % F in der characteristischen gleichung überein. Die zweite Lösung ist
     % nicht relevant. Wurde für mehrere Fälle überprüft define the
     % characteristic equation
     
% For symmetric and Hermitian matrices, the eigenvalues and singular values
% are obviously closely related. A nonnegative eigenvalue,  ??0, is also a singular
% value,? = ? . The corresponding vectors are equal to each other, u = v = x. A
% negative eigenvalue, ? < 0, must reverse its sign to become a singular value,
% ? = | ? | . One of the corresponding singular vectors is the negative of the other,
% u = ? v = x 

     cf = sym( det( F - eye(3) ) ) == 0;
     
     solutions = double( solve( cf, g) );   %, 'Real', true);

     % If for this equation no solution could be found assign 8
     if isempty( solutions )
         shear_amount( nr,1) = 8;
         shear_amount( nr,2) = i;
         nr = nr +1;
         continue
     end
     % if the solution for g is very high or low respectively set it to 7
     % or -7 for a good readability of the results
     %     for si = 1:2
     solution = solutions;
     if solution > 10
         shear_amount(nr,1) = 7;
     elseif solution < 10
         if solution < -10
             shear_amount(nr,1) = -7;
         elseif solution > - 10
             shear_amount(nr,1) = solution;
         end
     end
     shear_amount(nr,2) = nr;
     
     %% Transform the modified Bain strain matrix back into the cubic system
     % for a direct comparison between all martensite variants
     
     G1 = [1 0 shear_amount(i); 0 1 0; 0 0 1];
     U_j(:,:,nr) = R' * ( U_rot_m_l * G1 ) * R;
     eigenvalues(nr,:) = eigs( U_j(:,:,nr) );
     [W,D] = eig( U_j(:,:,nr) );
     WW(:,:,nr) = W;
     nr = nr +1;
     %end
     
     %% if more than one solution for g is available
     %if ( l(i,:) == [1 -1 1] & m(i,:) == [-1 0 1] ) | ( l(i,:) == [1 1 1] & m(i,:) == [1 0 -1] )
         %U_j(:,:,nr-1)    
     %end
     
     % if length(gg) == 2
     %     g2 = double( gg(2) )
     %     G2 = [1 0 g2; 0 1 0; 0 0 1];
     %     B_mod2 = G2 * U_rot;
     %     eigs( B_mod2 )
     % end

 end
 
 % short shear amout for ascending magnitude of shear
[~,idx] = sort( abs(shear_amount),1 )
shear_amount = shear_amount( idx ) 
 
%% Write data of minimum values for g to a file
format_variant = ' U%d = %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n \n';
 fprintf(fileID,format_variant,variantnr,B);
 
 format_m_aust = 'm_aust   = (%d, %d, %d) \t';
 format_l_aust = 'l_aust   = [%d, %d, %d] \t';
 format_m_mart = 'm_mart = (%d, %d, %d) \t';
 format_l_mart = 'l_mart = [%d, %d, %d] \n';
 
 format_g = 'g = %7.4f \t';
 format_e = 'eigs = %7.4f %7.4f %7.4f \n';
 format_mod_variant = ' U_mod = %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n \n';

 mJa = aJm \ eye(3);
for i = 1:size(eigenvalues,1)
% in case the shear is formally applied in the parent phase
    m_aust = m( idx(i,1), : );
    l_aust = l( idx(i,1), : );
    l_mart = mJa*l_aust';
    m_mart = m_aust * aJm;

% in case the shear is formally applied in the product phase
%     m_mart = m( idx(i,1), : );
%     l_mart = l( idx(i,1), : );
%     l_aust = aJm*l_mart';
%     m_aust = m_mart * mJa;

    g = shear_amount(i);
    e = eigenvalues( idx(i,1),:);
    mod_variant = U_j(:,:,idx(i,1) );
    % ev = eigenvectors =
    
    fprintf(fileID,'shear system nr: %d \n', idx(i,1) );
    fprintf(fileID,format_m_aust, m_aust);
    fprintf(fileID,format_l_aust, l_aust);
    fprintf(fileID,format_m_mart, m_mart);
    fprintf(fileID,format_l_mart, l_mart);
    fprintf(fileID,format_g, g);
    fprintf(fileID,format_e, e);    
    fprintf(fileID,format_mod_variant, mod_variant); 
end
fprintf(fileID,'\n ---------------------------------------------- \n \n');
fclose('all');
%%
min_g_U_j_out = min_g_U_j_in;
WW_out = WW_in;
for i = 1:size(U_j,3)
    min_g_U_j_out(:,:,(size(min_g_U_j_out,3)+1) ) = U_j(:,:,idx(i,1) );
    WW_out(:,:, size(WW_out,3)+1 ) = WW(:,:,idx(i,1) ); 
end
 

  
 end
