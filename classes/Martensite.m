classdef Martensite < Base 
    
    properties(Access = public)
        U; % = B - Bain strain / structural-/Transformation stretch tensor - hermitian part of the deformation gradient
        F; % Deformation_gradient of transformation - To determine hermitian part via polar decomposition e.g. for NiTi
        % e_mart = F * e_aust
        C_am; % lattice correspondence (similar to OR and comprising U!)
        % TODO - integrate function -  correspondence_matrix_components(lattice1, lattice2, m1, m2, use_planes)
        % to determine lattice correspondence - or see X. Chen Paper or
        % their Code - for now we take the hard coded Bain-correspondence from the well
        % known figure...
        %R % rotational part of the
        IPS_solutions; % = Solution_array();
        ILS_solutions;
        block_solutions;
        %
        considered_plasticity; % 1-only mart slip systems, 2-only austenite slip systems, 3-slip systems of both lattices
        %
        invariant_lines;
    end
    properties (Dependent)
        cp;
    end
    
    methods
        % constructor
        function obj = Martensite()
            obj.IPS_solutions = Solution_array();
            obj.ILS_solutions = Solution_array();
            %obj.block_solutions = Solution_array_block
        end
        %------------------------------------------------------------------
        %TODO implement loadobject
        % function obj = joadobj( file )
        %------------------------------------------------------------------
        function obj = set.U(obj, B_in)
            if abs(B_in - B_in') < 1.e-9
                if eig( B_in ) > 0
                    obj.U = B_in;
                else
                    error('Specified Bain strain not positive definite')
                end
            else
                error('Specified Bain strain not symmetric')
            end
        end
        %-------------------------------------------------------------------
        function U = get.U(obj)
            % if the matrix is symmetric and positive definite
            % (A matrix is positive definite if all its associated
            % eigenvalues are positive)
            if ( (sum(sum(abs(obj.U - obj.U'))) < 1.e-9 ) && (all(eig( obj.U )) > 0.) ) % if U is already set
                U = obj.U;
            else
                U = polardecomposition(obj.F);
            end
        end
        %-------------------------------------------------------------------
        function obj = set.considered_plasticity(obj, nr)
            if ismember (nr, [1,2,3])
                obj.considered_plasticity = nr;
            else
                error('martensite.considered_plasticity has not been set properly (can only be 1, 2 or 3) see Martensite class')
            end
        end
        %-------------------------------------------------------------------

%  09.05.2017----TODO - remove or adopt so that it fits with IPS_solution
% 
%        function obj = set.F(obj, F)
%             % In the case of martensitic transformations, a further condition has to be satisfied;
%             % the lattice transformation strain must also be an invariant line strain if the interface is
%             % to be glissile (see christian, crocker - dislocations and lattice transformations - Dislocations in crystals vol 3).
%             % Also F cannot be uniquely determined because there is an infinite number of
%             % Lattice correspondances ai = F*bi.
%             % Moreover, we limit ourselves to deformations that preserve orientation,
%             % i.e. those with det(F) > 0 (tripe products have the same sign)
%            if det(F_in) < 0.
%                error('det(F) < 0, non-orientation preserving transformation')
%            else
%                obj.F = F;
%            end
%        end
%%-------------------------------------------------------------------
%         function def = get.F(obj)
%             display('testtest')
%             if abs(obj.F - eye(3)) < 1.e-5
%                 def = eye(3);
%             else
%                 def = obj.B * obj.R;
%             end
%         end
%%------------------------------------------------------------------
%         function def = F_from_atom_postions( lattice1, lattice2, A1, A2 )
%          Programm if needed, maybe not very reasonable...     
%         end
%%------------------------------------------------------------------     
%         function getvariant( index ):
%         % given the variant index return the according variant
% 
%         if self.__Ulist is None:
%             self.calcvariants()
%         return self.__Ulist[n - 1]  
%------------------------------------------------------------------
        function cp = get.cp(obj)
            cp = obj.U * obj.C_am;
        end
%-------------------------------------------------------------------
        function vars = symmetry_related(obj, variant)
            % calculate all symmetry related variants of the transformation
            % from one initial "variant"
            % store the matrices associated to variant_i
            % and the indices of elements in the Laue group R_i that maps
            % the initially given U1 to them. Also store Rotations that
            % collapse to the same variant (if any)
            if sum(size( variant )) == 4 % size(variant) = [1,3]
                for i=1:size(obj.Point_group.matrices,3)
                    vars(:,i) = obj.Point_group.matrices(:,:,i) * variant; % active rotation
                end
            elseif sum(size( variant )) == 12; %[3,3,3,3] % For 4th order tensor
                rotateTensor4(obj.Point_group.matrices(:,:,i), variant)
            else
                for i=1:size(obj.Point_group.matrices,3)
                    vars(:,:,i) = obj.Point_group.matrices(:,:,i) * variant * obj.Point_group.matrices(:,:,i)';
                end
            end
        end
        %-------------------------------------------------------------------   
        %         TODO rewrite to matlab - structrans (Phyton Code)
        function bool = isreversible( obj )
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
    end % methods end
end % class end

%         function obj = set.R(obj, R_in)
%             %R_in*R_in'
%             %det(R_in)
%             if in_O3(R_in);
%                 obj.R = R_in;
%             else
%                 error('Specified Rotation is not a proper orthogonal Matrix')
%             end
%         end
        %------------------------------------------------------------------
%         function rot = get.R(obj)
%             if abs(obj.F - obj.F') < 1.e-9  && eig( obj.F ) > 0
%                 [~, rot] = polardecomposition(obj.F);
%             end
%         end
        %------------------------------------------------------------------
        

   
% Verify - The following probably only holds if the basis remains the same:
% function bool = unextended( vec, B )
% % given a vector and the Bain strain in the same Basis, this function
% % verifies wheter the vector remains undistorted by the Bain strain
% see Bhadeshia- Worked examples in the Geometry of crystals p.8
% if abs( (B(1,1)^2 -1.)*vec(1) + (B(2,2)^2 -1.)*vec(2) + (B(3,3)^2 -1.)*vec(1) ) < 1.e-9
%     bool = true
% else
%     bool = false
% end
