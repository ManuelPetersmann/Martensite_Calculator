% Everything due to WOLFGANG PRANGER Masterthesis, same coord syst as Waitz - Size
% effects paper (who uses lattice parameters as in Knowles 1981)
% Lattice constants due to Kudoh et al.
parent =  Martensite();
parent.my_base = eye(3);
product = Base();
product.my_base = eye(3);
% for the product phase no Martensite object is necessary since all calculations
% are carried out in the coordinate system of the parent phase
product.Bravais_type = 'monoclinic';
parent.Bravais_type = 'cubic';
parent.Centering = 'P';
product.Centering = 'P';
product.Point_group; 

beta = degtorad(97.78);
ac = 3.015;
am = 2.898;
bm = 4.108;
cm = 4.646;
parent.Lp(1:3) = ac; %{1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
           % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
product.Lp(1:6) = [am, bm, cm, pi/2, beta, pi/2];

% components martensite
m1 = [1 0 0;
      0 1 0;
      0 0 1];
%components austenite
m2 = [1  0 0;
      0  1 1;
      0 -1 1];
or = lattice_correspondance();
cp = or.correspondance_matrix_components(product, parent, m1, m2)
 


% cartesian system of the monoclinc cell = basis of the tetragonal cell
Fm = [ am/ac            0          cm/(sqrt(2)*ac) * cos(beta);
      0         bm/(sqrt(2)*ac)       0;
      0                0          cm/(sqrt(2)*ac) * sin(beta)];
  
% Assemble TransformationMatrix for the monoclinic to cubic transformation
% Hane and shield % negative rotation around x-axis
R_to_cubic = [ 1 0 0; 0 cos(pi/4) -sin(pi/4); 0 sin(pi/4) cos(pi/4) ]';

% transformation to cubic system ( same as for tensor from monoclinic with
Fc = R_to_cubic * Fm * R_to_cubic';

%parent.F = Fc;

% koordinate transformation from CS used in my DA to that by Waitz (See
% hane and shield )
% if to_waitz == 1
%     Q =  [0 1 0; 0 0 1; 1 0 0]';
% else
%     Q =  [0 1 0; 0 0 1; 1 0 0];
% end

