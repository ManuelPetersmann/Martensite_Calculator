function [ y1, y2, y3, e1, e2, e3] = sorted_eig_vals_and_vecs( Ct )
% This function takes a matrix and returns its eigenvalues in ascending 
% order as well as the corresponding (normalized!) eigenvectors 

[V,D] = eig( Ct );
% D... Diagonal matrix of eigenvalues
% V... colum vectors are correspondig right eigenvectors i.e. A*V = V*D
eigs = [D(1,1),D(2,2),D(3,3)];
[eigs_sort, old_idx] = sort(eigs);
% arrange eigenvectors to order of sorted eigenvectors
y1 = eigs_sort(1);
y2 = eigs_sort(2);
y3 = eigs_sort(3);
e1 = V(:,old_idx(1));
e2 = V(:,old_idx(2));
e3 = V(:,old_idx(3));
e1 = e1 / norm(e1);
e2 = e2 / norm(e2);
e3 = e3 / norm(e3);

end

