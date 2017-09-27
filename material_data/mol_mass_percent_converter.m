function conv_vec = mol_mass_percent_converter( fracs, stringInput, to_molpercent )
% This function takes a 2xN cell-array 'fracs'. The first line holds
% weight-% or atom-% (mol fractions) (between zero and one) and the second 
% according molmasses or strings of the elements e.g. 
% fracs = {0.1, 16.0;   9.0,  58.69} or {0.1, 16.0; 'C','Ni'}
% returns a vector of the other fraction e.g.   [xi_C   xi_Ni] -> [w_C  w_Ni]
% stringInput: a boolean controlling determining wheter the
% Molemasses are input directly or as a Vector of Strings of the Element
% labels according to the periodic table
% to_molepercent - determines wheter mole or weightpercent are calculated

% check if all fractions sum to 1 with 1% tolerance
total = 0.;
for j= 1:size(fracs,2)
    total = total + fracs{1,j};
end
if abs(total - 1) > 0.01
    error('fractions do not sum to one please fix that.')
end

%% Get Molemasses if second line of fracs is given as Element-strings
if stringInput
    for j= 1:size(fracs,2)
        fracs{2,j} = molMass( fracs{2,j} );
    end
end

%%
denominator = 0.;
if to_molpercent
    % conv = xi ... mol percent
    for j= 1:size(fracs,2)
        conv = fracs{1,j} / fracs{2,j};
        if denominator == 0.
            for i = 1:size(fracs,2)
                denominator = denominator + fracs{1,i} / fracs{2,i};
            end
        end
        conv_vec(j) = conv / denominator;
    end
else
    % conv = wi ... mass percent
    for j= 1:size(fracs,2)
        conv = fracs{1,j} * fracs{2,j};
        if denominator == 0.
            for i = 1:size(fracs,2)
                denominator = denominator + fracs{1,i} * fracs{2,i};
            end
        end
        conv_vec(j) = conv / denominator;
    end
end


end

