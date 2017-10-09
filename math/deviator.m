function [ dev ] = deviator( A )
% call: deviator( A )
% Given a second order tensor in the form of a matrix
% calculates its deviator

if size(A,1) ~= size(A,2)
    error('matrix is not symmetric')
end

dev = A;
trA = (1./3.)*trace(A);
for i = 1:size(A,1)
        dev(i,i) = dev(i,i) - trA;
end


end

