function M_ave = average_molar_mass( fracs, stringInput )
% fracs {0.3 0.7;
%        55.4, 32.8}

if stringInput
    for j= 1:size(fracs,2)
        fracs{2,j} = molMass( fracs{2,j} );
    end
end

M_ave = 0.;

for i = 1:size(fracs,2)
    M_ave = M_ave + fracs{1,i} * fracs{2,i};
end

end

