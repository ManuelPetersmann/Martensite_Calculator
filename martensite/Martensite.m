classdef Martensite < Lattice
    
    properties
        F % Deformation_gradient of transformation
    end
    properties (Dependent)
        B % Bain_strain
    end
    
    methods
        % constructor
        function obj = Martensite()
        end
        %-------------------------------------------------------------------
        %TODO implement loadobject
        % function obj = joadobj( file )
        %-------------------------------------------------------------------
        %function obj = set.F(obj, F_in)
        %    
        %end
        %-------------------------------------------------------------------
        function U = get.B(obj)
            if is_symmetric(obj.F) && is_positive_definite(obj.F)
                U = obj.F;
            else
                U = polardecomposition(obj.F);
            end
        end
        %-------------------------------------------------------------------
        function getvariant( index ):
        % given the variant index return the according variant

        if self.__Ulist is None:
            self.calcvariants()
        return self.__Ulist[n - 1]        
        %-------------------------------------------------------------------   
        function [] = calcvariants( U, Point_group)
        % calculate all the variants of the transformation.
        % store the matrices associated to each variant U_i
        % and the indices of elements in the Laue group R_i that maps
        % the initially given U1 to them. Also store Rotations that
        % collapse to the same variant (if any)

        if self.__Ulist is None:
            if not self.isreversible():
                raise AttributeError(
                    "variants for irreversible martensitic transformations have not been implemented yet"
                )

            U1 = self.getU()
            ulist = [U1]
            idx = [[]]
            for i, Q in enumerate(self.__Laue.matrices()):
                V = Q.dot(U1).dot(Q.T)
                newU = True
                for j, U in enumerate(ulist):
                    if np.allclose(U, V):
                        newU = False
                        idx[j].append(i)

                if newU:
                    ulist.append(V)
                    idx.append([i])

            self.__Ulist = np.array(ulist)
            self.__Laueidx = np.array(idx)
        %-------------------------------------------------------------------   
        
            
        %TODO        
        %         TODO rewrite to matlab
        %            def isreversible(self)
        %
        %                 It is reversible if the symmetry group of ``U``
        %                 is a proper subgroup of the object's ``Laue`` group
        %
        %                 return whether the transformation is reversible
        %                 rtype boolean
        %                 raises AttributeError
        %
        %                         if `U` has not been assigned
        %
        %                 if self.__U is None
        %                     raise AttributeError(U has not been initialized)
        %                 lg = Lattice(self.getU()).getLauegroup()
        %                 lg0 = self.getLaue()
        %                 return lg.order() = lg0.order() and lg0.hassubgroup(lg)
        %                 if self.__U is None
        %                     raise AttributeError(U has not been initialized)
        %                 lg = Lattice(self.getU()).getLauegroup()
        %                 lg0 = self.getLaue()
        %                 return lg.order() = lg0.order() and lg0.hassubgroup(lg)
        
    end
end