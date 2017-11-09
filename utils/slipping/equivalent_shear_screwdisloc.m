function [screw_syst] = equivalent_shear_screwdisloc(s,m)
% call: equivalent_shear_screwdisloc(s,m)
% s || b ... shear direction, burgers vector of edge dislocation -> normal
% to line segment of dislocation for edge disloc
% return burgersvector of screwdisloc - b \perp s - that gives the same continuum shear

b = cross(s,m);

screw_syst = [mat2str(b) ,' ', m];



end

