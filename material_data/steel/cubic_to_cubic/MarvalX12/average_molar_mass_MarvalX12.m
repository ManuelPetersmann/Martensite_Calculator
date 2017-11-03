fracs = {0.1215, 0.0905, 0.0203, 0.007, 0.0035, 0.0;
         'Cr',   'Ni',   'Mo',   'Al',  'Ti',   'Fe'};
fracs{1,6} = 1. - sum(cell2mat(fracs(1,1:5))) % Calculate Fe fraction

average_molar_mass( fracs, 1 )