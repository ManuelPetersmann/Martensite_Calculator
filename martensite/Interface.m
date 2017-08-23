classdef Interface
    
    properties
        v % probe vector of interface
    end
    
    methods
        function Bc = Frank_Bilby(S1, S2, v, ref_to_new)
            % dislocation content in interface between deformations S1 and S2
            % depending on the definition of reference in the direction of a 
            % probe vector v [3,1] lying in the interface
            if ref_to_new
                Bc = (inverse(S1) - inverse(S2) ) * v;
            else
                Bc = (S1 - S2) * v;
            end
        end
    end
    
end

