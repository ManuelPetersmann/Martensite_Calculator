function[bool] = equal_matrices(mat1, mat2, epsilon)
size1 = size(mat1);
size2 = size(mat2);
count = 0;
if size1(1) == size2(1) && size1(2) == size2(2)
    for i = 1:size1(1)
        for j = 1:size1(2)
            if abs( mat1(i,j) - mat2(i,j) ) < epsilpon
                count = count +1;
            end
        end
    end
    if count == size1(1)*size1(2)
        bool=1; %entspricht true
    else
        bool=0; %entspricht false
    end
else 
    display('input matrizes must be the same size')
end

