function conv_vec = mol_mass_percent_converter( fracs, stringInput, to_molpercent )
% This function takes an array fracs composed of two vectors, the first holds
% weight% or atom% (mol fractions) (between zero and one) and the second 
% according molmasses or strings of the elements e.g. 
% fracs = {0.1, 16.0;   9.0,  58.69} or {0.1, 16.0; 'Fe','Ni'}
% returns a vector of the other fraction e.g.   [xi_C   xi_Ni] -> [w_c  w_Ni]
% frac = Nx2
% stringInput, to_molepercent are booleans controlling determining wheter the
% Molemasses are input directly or as a Vector of Strings of the Element
% labels according to the periodic table

if stringInput
    for j= 1:size(fracs,1)
        fracs{2,j} = molMass( fracs{2,j} );
    end
end

if to_molpercent
    % conv = xi ... mol percent
    for j= 1:size(fracs,2)
        conv = fracs{1,j} / fracs{2,j};
        for i = 1:size(fracs,2)
            denominator = denominator + fracs{1,i} / fracs{2,i};
        end
        conv_vec(j) = conv / denominator;
    end
else
    % conv = wi ... mass percent
    for j= 1:size(fracs,2)
        conv = fracs{1,j} * fracs{2,j};
        for i = 1:size(fracs,2)
            denominator = denominator + fracs{1,i} * fracs{2,i};
        end
        conv_vec(j) = conv / denominator;
    end
end


end

