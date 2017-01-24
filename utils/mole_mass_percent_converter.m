function [ conv ] = mole_mass_percent_converter( frac_molemass, to_molpercent )
% This function takes an array frac_molemass with weight or mole fractions
% (between zero and one) and according Molemasses e.g. frac_molemass = [[0.1, 16.0], [9.0,   58,69] ... ] 
% of each element and returns the other fraction                        xi_C     M_c   xi_Ni   M_Ni
% of the first element in the frac_molemass array respectively.
% dimension frac molemass = 1xNx2


if to_molepercent
    % conv = xi ... mole percent
    conv = frac_molemass(1,1)/frac_molemass(1,2);
    for i = 1:size(frac_molemass,1)
        denominator = denominator + frac_molemass(i,1)/frac_molemass(i,2);
    end
    conv = conv / denominator;
else
    % conv = wi ... mass percent
    conv = frac_molemass(1,1) * frac_molemass(1,2);
    for i = 1:size(frac_molemass,1)
        denominator = denominator + frac_molemass(i,1) * frac_molemass(i,2);
    end
    conv = conv / denominator;
end


end

