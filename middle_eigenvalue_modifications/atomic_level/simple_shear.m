function [solutions] = simple_shear(B, cp, ns, ds )

 for i = 1:size(ms,1) % since there are two solutions for g for each slip system respectively
     d(i,:);
     n(i,:);

     %% calculate the representation of the Bain strain Tensor in the shear
     % coordinate system. Therefore get the transformation matrix that
     % transforms the specified U_i into the basis of the shear system:
     R = Basiswechselmatrix( [1 0 0], [0 1 0], ( d(i,:)/ norm( d(i,:) ) ), ( n(i,:)/ norm( n(i,:) ) ) );
     U_rot_m_l = R * B * R'; % 
     
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
     F = sym( U_rot_m_l * G );     

     cf = sym( det( F - eye(3) ) ) == 0;
     
     solutions = double( solve( cf, g) );   %, 'Real', true);

     % If for this equation no solution could be found assign 8
     if isempty( solutions )
         continue
     end
     % if the solution for g is very high or low respectively 
     if (solution > g_max) || (solution < (-1)*g_max)
         continue
     end
     
     %% Transform the modified Bain strain matrix back into the cubic system
     % for a direct comparison between all martensite variants
     
     G1 = [1 0 shear_amount(i); 0 1 0; 0 0 1];
     U_j(:,:,nr) = R' * ( U_rot_m_l * G1 ) * R;
     eigenvalues(nr,:) = eigs( U_j(:,:,nr) );
     [W,D] = eig( U_j(:,:,nr) ); % eigenvalues and eigenvectors
     WW(:,:,nr) = W;
     nr = nr +1;
     % end

 end
 
 % short shear amout for ascending magnitude of shear
[~,idx] = sort( abs(shear_amount),1 )
shear_amount = shear_amount( idx ) 
 
  
 end
