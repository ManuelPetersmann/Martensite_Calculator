classdef Base < Bravais_Lattice
    % This class provides the user defined lattice base and transformation functions
    % for coordinate transforms. 
    % A Basis is stored in a 3x3 matrix as colum vectors
    % Vectors must be given as colum vectors
    % Note that the baseis (=bases) here are generally not cartesian but affin!
    % The words components of a vector and coordinates of a point are used synonimous
    % International tables of crystallography Vol A, A1
    properties (SetAccess = public)
        my_base; % can e.g. be a cartesian (orthonormal) coordinate system so that everything
        % calculated in it can be directly used as input in e.g. finite element calculations
        % without specifying a basis change.
    end
    properties (Dependent)
        Cov_metric;
        Contra_metric;
        Dual_base; % reciprocal = covariant basis = dual basis
        % for cubic crystal: components = miller indices of a plane whose normal is along that
        % vector, with the spacing of the plane given by the inverse of the magnitude of that vector.
    end
    %----------------------------------------------------------------------
    
    methods   
        % constructor
        function obj = Base( ) %basis )
            % obj.my_base = basis;
        end
        %------------------------------------------------------------------
        % function signatures
        metric = calc_metric(basis);
        db = calc_dual_basis(basis);
        %------------------------------------------------------------------
        % value class set functions must return the modified object to the
        % calling function. Handle classes do not need to return the
        % modificed object
        function obj = set.my_base(obj, components)     
            if size(components,1) ~= size(components,2)
                error('Conventional basis must be a square matrix');
            elseif det(components) <= 0.0
                error('Base matrix is singular!')
            else
                obj.my_base = components;
            end
        end
        %-------------------------------------------------------------------
        function metric = get.Cov_metric( obj )
            % given a covariant or a contravariant basis returns its metric
            metric = calc_metric(obj.my_base);
        end
        %-------------------------------------------------------------------
        function cont_metric = get.Contra_metric( obj )
            % given the basis vectors E of the primitive crystallographic lattice
            % calculate the corresponding Contravariant metric Z^i
            % mf = matrix_functions;
            cont_metric = inv(obj.Cov_metric);
        end
        %------------------------------------------------------------------
        function db = get.Dual_base( obj )
            % given a basis (co- or contravariant) returns its dual basis
            db = calc_dual_basis( obj.C );
        end
        %------------------------------------------------------------------                
        function bool = generate_same_lattice( B1, B2 ) % volume preserving but including shears
            % given two equally oriented lattice bases (det > 0 -> right handed systems)
            % as colum vectors in two matrices, this function calculates 
            % wheter they generate the same lattice % see e.g. 
            % Bhattacharya - Microstructures of Martensites p.32
            % Any two sets of primitive lattice vectors for a given lattice are related by
            % a lattice invariant transformation i.e., a unimodular matrix
            % of integers. 
            %matfuncs = matrix_functions();
            mu = B2 * inverse(B1);
            if abs(det(mu) - 1.) <  1.e-9
                bool = true;
            else
                bool = false;
            end
        end
%            def lattice_group(self):
%             E = self.__E % conventional basis
%             Einv = la.inv(E)
%             mats = np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.point_group().matrices()], dtype='int') 
%             return MatrixGroup(mats)
        %-------------------------------------------------------------------        
        function bool = is_isometry( P, is_reciprocal )
            % Given an affine Mapping M and a basis (real -> second argument = true
            % or reciprocal -> second argument = false) calculates wheter for the given
            % basis the mapping is an isometry (distances and angles are
            % unchanged i.e. the metric remains the same).            
            % Note: Two sets of lattice vectors that are related by a
            % rotation generate the same metric
            % matfuncs = matrix_functions;
            Q = inv(P);
            if is_reciprocal
                if abs(Contra_metric - Q * Contra_metric * Q') < 1.e-9 % see e.g. Int. table cryst vol A, p.85
                    bool = true;
                else
                    bool = false;
                end
            elseif abs(Cov_metric - P' * Cov_metric * P) < 1.e-9 
                bool = true;
            else
                bool = false;
            end
        end
        % TODO check if above and below functions are equivalent
        function in_symmetry_group( obj, H )
            % the symmetry- or lattice-group is the set of deformations that map a lattice back to itself
            % this group also contains shears. The pointgroup therefore is a subgroup of this group
            % %see Bhattacharya MM - p.33
            %             for i=1:size(lattice.)H
            %             end
        end        
        %------------------------------------------------------------------
    end % end methods
end % end class

