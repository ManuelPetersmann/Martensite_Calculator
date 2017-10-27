function [det_eq, grad_det] = det_grad_sym(F1, F2, x, U)
%
syms FC(y)
FC(y) = y * F1  +  (1.-y) * F2

det_eq = det( FC(y) )

grad_det = diff(det_eq,y)

grad_det = subs(grad_det,x)

end

