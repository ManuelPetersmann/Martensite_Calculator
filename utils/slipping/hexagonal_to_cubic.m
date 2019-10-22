function [uc, vc, wc] = hexagonal_to_cubic_indices( uh, vh, wh, th)
% call: hexagonal_to_cubic_miller(U, V, W, T)

uc = 2*uh + vh;
vc = uh + 2*vh;
wc = wh;

end

