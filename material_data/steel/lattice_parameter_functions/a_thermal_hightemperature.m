function [a_austenite, a_ferrite] = a_thermal_hightemperature( T )
%A_GAMMA_THERMAL 
% J.H. Root, M. Onink, C.M. Brakman, F.D. Tichelaar, S. Van der Zwaag, E.J.
% Mittemeijer, and N.B. Konyer. The lattice-parameters of austenite and ferrite
% in Fe-C alloys as functions of carbon concentration and temperature. Scripta
% Metallurgica Et Materialia, 29(8):1011{1016, 1993.

% temperature intervall 1030 - 1250 K
a_austenite = 0.36511 (1 + 23.3e-6*(T-1000))
% ferrite lattice parameters 800 K to approximately 1050 K 
a_ferrite = 0.28863*(1 + 17.5e-6*(T - 800));
end

