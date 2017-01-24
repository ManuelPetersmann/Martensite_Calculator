function epsilon_c = constrained_strain_eshelby( S, epsilon_t )
% call: constrained_strain_eshelby( S, epsilon_t )
% S_ijkl ... Eshelby tensor
% epsilon_t ... transformation strain

epsilon_c = zeros(3);

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                epsilon_c = epsilon_c + S(i,j,k,l) * epsilon_t(k,l);
            end
        end
    end
end

end

