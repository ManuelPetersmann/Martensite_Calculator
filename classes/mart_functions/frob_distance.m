function d = frob_distance(A,B)
%
C = A - B;
d = sqrt(trace( C*C' ) );

end

