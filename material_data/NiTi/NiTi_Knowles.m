% EVERYTHING DUE TO Knowles 1981
ac = 3.015
cm = 4.646 %am me
am = 2.898 %bm me
bm = 4.108 %cm me
beta = degtorad(97.78)
F = (1/ac) * [bm/sqrt(2)         0                  0;
                0           cm*sin(beta)/sqrt(2)    0;
                0           cm*cos(beta)/sqrt(2)    am]
            
% Knolws 1981, rotation around z-axis
R_to_cubic = assembleTransformationMatrix( [1 -1 0], [1 1 0], [1 0 0], [0 1 0] )
