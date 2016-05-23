function bool = in_vec_array( vec, vec_array)
% this function takes a 1xN vector - vec and an PxN matrix - vec_array
% and checks whether vec is contained in  vec_array

bool = false;
for j = 1:size(vec_array,1)
    if abs(vec - vec_array(j,:) ) < 1.e-4
        bool = true;
        break
    end
end

end

