function [] = plot_quadratic_surface(coefs, xyrange, points)

if length(coefs) ~= 10
    error('coef vector must have 10 numbers (coefficients of most general quadratic surface)');
end
if nargin < 2
    xrange = 50.;
    yrange = 50.;
end

coefCell = num2cell( coefs );
[A,B,C,D,E,F,G,H,I,J] = coefCell{:};

gv = linspace(-xyrange(1), xyrange(2), points); % adjust for appropriate domain
[xx, yy, zz]=meshgrid(gv, gv, gv);  
F = A* xx.^2 + B* 2*yy.^2 + C * zz.^2 + D*xx.*yy + E*xx.*zz + F*yy.*zz + G*xx + H*yy + I*zz + J;

figure
isosurface(xx, yy, zz, F, 0)

end

