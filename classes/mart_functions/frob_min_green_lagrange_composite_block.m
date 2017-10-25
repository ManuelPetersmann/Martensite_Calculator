function [x, d_old, Fcomp] = frob_min_green_lagrange_composite_block( F1, F2, x_ini) %, min_sener_change )
% call: x = frob_min_green_lagrange_composite( in )
% optimize |Fc^t * Fc - I| = frob(Green_lagrange) , where Fc is a linear mixture
% Fc =x1 F1 + x2 F2;   x1 +x2 = 1;   Fi are Matrices.
% using intervall sectioning
% Minimization proposed by Muehlemann2015 - Morphology of lath martensite a new perspective
% see Fig. 5 there.

min_sener_change = 1.e-5; % Has been found to be a reasonable value!
F1  = reshape(F1,9,1);
F2  = reshape(F2,9,1);
Fcoefs  = cat(2,F1,F2);
%
x1 = x_ini;
dx = 0.01;
%% determine search direction of start or do nothing if already at minimum
d_old = strain_ener_measure(x1); 
d_new1 = strain_ener_measure(x1+dx);
d_new2 = strain_ener_measure(x1-dx);
if ((d_old - d_new1) > min_sener_change) || ((d_old - d_new2) > min_sener_change)
    if d_new1 > d_new2
        dx = -1*dx;
    end
else
    x = [x1; 1.-x1];
%    d = d_old;
    Fcomp = Fcoefs*x;
    return
end

%iters = 1;
%% continue to search in initial direction until convergence
while true
    d_new = strain_ener_measure( x1+dx );
    if d_new > d_old % if residual increases cut back
        dx = 0.5*dx;
    else
        x1 = x1 + dx; % if residual decreases update x1
        d_old = d_new;
    end
    if abs(d_old - d_new) < min_sener_change % if increment change falls below a threshold - stop 
        break   
    end   

%    iters = iters + 1;
end
%iters
x = [x1; 1.-x1];
Fcomp = reshape(Fcoefs*x ,3,3);

%c = linspace(0,1,1000);
%d = arrayfun(@strain_ener_measure, c);
%plot(c,d)

    function d = strain_ener_measure(x)
        y = 1. - x;
        xi = [x; y];
        Fc = Fcoefs*xi;
        Fc = reshape(Fc,3,3);
        E = Fc'*Fc - eye(3);
        d = sum(dot(E,E));
    end

end