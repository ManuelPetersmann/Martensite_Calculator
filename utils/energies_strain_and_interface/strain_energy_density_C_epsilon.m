function SENER = strain_energy_density_C_epsilon( C, epsilon )
% call: strain_energy_density_C_strain( C, epsilon )
% calculates the strain energy density of an anisotropic, linear elastic
% solid given the elastic constant tensor C (3x3x3x3) and the strain tensor
% epsilon (3x3)
% U = 1/2 * C_ijkl * epsilon_ij * epsilon_kl
% U = 1/2 * S_ijkl * sigma_ij * sigma_kl

% strain_energy_density - SENER = C_ijkl epsilon_ij epsilon_kl 

SENER = 0;


for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                SENER = SENER + C(i,j,k,l) * epsilon(i,j) * epsilon(k,l);
            end
        end
    end
end
SENER = SENER *0.5;


end

