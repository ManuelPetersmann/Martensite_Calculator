function P = basis_relation( basis1, basis2 )

P = basis1 * inv( basis2 );
% Difference between here and Ulrich Müller:
% He writes: basis1 = basis2 * P -> the Basis Vectors are assembled as row
% vectors in the matrices. Here it is written: basis1=P*basis2 -> the Basis
% Vectors are assembled as colum vectors in the matrices. 
% Note (AB)' = B'A'
% Also written as [B1 |> B2] see Cayron 2006 - Groupoid of orientational
% variants
end