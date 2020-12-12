ac = 3.015
am = 4.646 % equal to cm knowles
bm = 2.898 % equal to bm knowles
cm = 4.108 % equal to am knowles
gamma = degtorad(97.78);

F = [ sin(gamma)*am/(sqrt(2)*ac)    0       0;
      cos(gamma)*am/(sqrt(2)*ac)    bm/ac   0;
      0             0         cm/(sqrt(2)*ac) ]
  
% cartesian system of the monoclinc cell = basis of the tetragonal cell
Xm = [1, 0, 0];

Ym = [0, 1, 0];

Zm = [0, 0, 1];


% basis of the cubic cell 
Xc = [1/sqrt(2), 0, 1/sqrt(2)];

Yc = [0, 1, 0];

Zc = [- 1/sqrt(2), 0, 1/sqrt(2)];


% Assemble TransformationMatrix for the monoclinic to cubic transformation
R_to_cubic = assembleTransformationMatrix( Xm, Ym, Xc, Yc );

% transformation to cubic system ( same as for tensor from monoclinic with
Fc = R_to_cubic * F * R_to_cubic';

