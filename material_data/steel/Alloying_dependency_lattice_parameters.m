function [ c_mart, a_mart, a_aust ] = Kurdjumov_Kaminsky_relations( a0_aust, cc )
%KURDJUMOV_KAMINSKY_RELATIONS from Turteltaub, Suiker 2005 JMPS p.16
% Lattice parameters of martensite as a function of carbon concentration (wt%)
% "cc" and of the lattice parameter of the austenite "a0_aust"
% Only valid for cc > 0.18 wt%

k1 = 0.116
k2 = 0.013
k3 = 0.044

a_aust = a0_aust + k3*cc;

a_mart = a0_aust - k2*cc;

c_mart = a_0m + k1*cc;

end

function [ c_over_a ] = Kurdjumov_Kaminsky_relations( cc )
% Honda - Nishiyama - martensite tetragonality relation
% cc - carbon concentration (wt%) - only valid for   cc > 0.6%
c_over_a = 1. + 0.0045;
end

function [ c_over_a ] = Kurdjumov_Kaminsky_relations( cc )
% tetragonality relation for 0-0.6 wt% carbon (cc -carbon content in wt%)
% after Lu Y., Richard - MSE-A 2017

c_over_a = 1. + 0.0031;
end



