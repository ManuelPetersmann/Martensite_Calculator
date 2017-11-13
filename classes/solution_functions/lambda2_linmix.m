function lambda2_mix = lambda2_linmix( x, F1, F2)
Fc = linmix2(x, F1, F2);
[~,lambda2_mix] = sorted_eig_vals_and_vecs(Fc'*Fc);
end