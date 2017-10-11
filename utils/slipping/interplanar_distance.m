function d = interplanar_distance( lp, pm, lattice )
% this function takes a vector lattice parameters lp= [a,b,c,alpha,beta,gamma] containing the lattice constants
% and angles depending on the Bravais lattice in this order but without
% non-used values, i.e. e.g. lp = [a,c] for tetragonal and 
% for monoclinic lp = [a,b,c,beta].
% The second argument is a vector of miller indices 'pm'(...plane miller) specifying the lattice
% plane for which the spacing should be calculated
% The third argument is a string identifying the Bravais lattice:
% Possible values are: 'cubic', 'tetragonal','hexagonal',rhombohedral',
% 'monoclinic','triclinic', 

switch lattice
    case 'cubic'
        x = sum(pm.^2) / lp(1);
    case 'tetragonal'
        x = (pm(1)^2 + pm(2)^2) / lp(1) + pm(3)^2 / lp(3);
    case 'hexagonal'
        x = 0.75*(pm(1)^2 + pm(1)*pm(2) + pm(2)^2)/lp(1) + pm(3)^2 / lp(3);
    case 'rhombohedral'
        x = sum(pm.^2)*sin(lp(4))^2 +  2*( pm(1)*pm(2) + pm(2)*pm(3) + pm(1)*pm(3) ) * ( cos(lp(4))^2 - cos(lp(4)) );
        x = x / ( lp(1)^2* ( 1. - 3*cos(lp(4))^2 + 2*cos(lp(4))^3 ) );
    case 'monoclinic'
        x = pm(1)^2/lp(1)^2  +  pm(2)^2*sin(lp(5))^2 / lp(2)^2  + pm(3)^2 / lp(3)^2 - 2*pm(1)*pm(3)*cos(lp(5)) / (lp(1)*lp(3));
        x = x*csc(lp(5))^2;
    case 'triclinic'
        % I guess I will never need it ^^ - https://en.wikipedia.org/wiki/Crystal_structure
end

d = 1/sqrt(x);

