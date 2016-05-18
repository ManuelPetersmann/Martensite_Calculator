classdef Bravais_Lattice
    % Every operation is in accordance with the International tables of
    % crystallography Vol A, A1
    
    properties (SetAccess = public)
        Bravais_type % a string specifying the Bravais lattice type e.g. 'cubic'
        Centering % a string specifying the exact type e.g. simple 'P',
        % face-centered 'F', body-centered 'I' etc.
        % Lattice_parameters
        Lp = [1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2]; 
    end
    properties (Dependent, Access = public)
        Point_group  % 3x3xN matrix containing the point group symmetry matrices
        % This group is often also denoted as Laue-Group, which however
        % also contains matrices generated in a left-handed system (det<0)
        % which we do not consider here
        E
        % primitive (crystallographic) basis vectors e_i stored as E(:,:,1) =
        % [e1;e2;e3], row vectors
        % i.e. lattice vectors with integer coefficients t_i: {t_vec = t1*a1 + t2*a2 + t3*a3}.
        % this means that the a_i are translation vectors!
        % This matrix array also contains the conversion matrix to the
        % conventional basis in E(:,:,2)
        C
        % conventional crystallographic = covariant = real basis = direct lattice
        % c_i stored as C = [c1;c2;c3], (such that as many angles as 
        % possible amout to 90°. e.g.
        % primitive lattice --> conventional lattice = primitive (e.g. primitive cubic),
        % otherwise "centered" e.g. face centered cubic.
        % note that for e.g. monoclinic, one conventional basis vector is
        % not parrallel to a primitive basis vector
        % These are the crystallographic bases used in the International Tables A !
        % must be user specified because of different possibilities
        % e.g. two for monoclinic
        Cov_metric_C
        Contra_metric_C
        reciprocal_base % = dual- or contravariant base
        Lattice_group
        density
    end
    
    methods
        % constructor
        function obj = Bravais_Lattice()
        end
        %------------------------------------------------------------------
        function obj = set.Bravais_type(obj, input)
            Bravais_strings = ['cubic       '; 'hexagonal   '; 'rhombohedral'; ...
                'tetragonal  '; 'tetragonal  '; 'orthorhombic'; 'monoclinic  '; ...
                'triclinic   '];
            Bravais_list = cellstr(Bravais_strings);
            if ismember( input, Bravais_list)
                obj.Bravais_type = input;
            else
                error('Invalid Bravais Type')
            end
        end
        %------------------------------------------------------------------
        function obj = set.Centering(obj, input )
            Centering_strings = ['P'; 'C'; 'I'; 'F'];
            if ismember( input, Centering_strings)
                obj.Centering = input;
            else
                error('Invalid Centering')
            end
        end
        %------------------------------------------------------------------
        %         TODO write function or class for user input e.g.
        %         bravais_type = input('choose bravais lattice') %later maybe make a
        %         GUI for that#
        function obj = set.Lp(obj, Lpars)
            %             switch obj.Bravais_type % the number of lattice parameters must match the bravais type
            %                 case 'cubic'
            %                     Par_num = [1, 0]; %number of different lengths and angles unequal 90°
            %                 case 'hexagonal'
            %                     Par_num = [2, 0];
            %                 case 'monoclinic'
            %                     Par_num = [3, 1];
            %                 case 'triclinic'
            %                     Par_num = [3, 3];
            %                     %TODO implement other lattices if needed
            %             end
            %             for i = 1 : obj.Par_num(2)
            %                 val = input( ['choose angle(s) because the', ...
            %                     'assignment of an coordinate system is not unique: alpha:', ...
            %                     'angle between a-b, beta b-c, gamma c-a ']);
            %                 if isnumeric(Lattice_params{i})
            %                     obj.Lattice_params{i} = val;
            %                 else
            %                     error('Value must be numeric')
            %                 end
            % deg2rad
            %             end
            obj.Lp = Lpars;
        end
        %------------------------------------------------------------------
        function point_group = get.Point_group( obj )
            point_group = matrix_group( obj.Bravais_type );
        end
        %------------------------------------------------------------------
        function bool = in_PointGroup(g)
            % Given another MatrixGroup g (3x3xN array), checks if its a subgoup of the
            % lattice point group
            % compare group order, i.e. number of elements in group
            if ~ismember( size(g,3), divisors( size(Point_group,3)) )
                bool = false;
                return
            else % check for matrices if they are contained
                for i=1:size(g,3)
                    if ~in_matrix_array( g(:,:,3), Point_group )
                        bool = false;
                        return
                    end
                end
                bool = true;
            end
        end
        %------------------------------------------------------------------
        function value = get.E(obj)
            % generate the base matrix for primitive cell of 3D Bravais lattices
            % it is stored in E(:,:,1)
            % in E(:,:,2) the conversion matrix W from the primitive P to the 
            % conventional base C is stored, i.e.:  C = W*P
            % See - e.g. Int. Table. Cryst. Vol A, p.81-83
            switch obj.Bravais_type
                case 'cubic'
                    switch obj.Centering
                        case 'P' % simple
                            value = obj.Lp(1) * eye(3);
                        case 'F' % Face centered fcc
                            value = 0.5 * obj.Lp(1)*[1, 1, 0;
                                                     0, 1, 1;
                                                     1, 0, 1]';
                            value(:,:,2) = 0.5* [1, 0, 1;
                                                 1, 1, 0;                  
                                                 0, 1, 1];
                        case 'I' % Body centered bcc
                            value = 0.5 * obj.Lp(1) * [1,  1, 1;
                                                      -1,  1, 1;
                                                      -1, -1, 1]';
                            value(:,:,2) = 0.5 * [1, -1, -1;
                                                  1,  1, -1;
                                                  1,  1,  1];
                    end
                case 'hexagonal'
                    value = cat2(obj.Lp(1) *[1,       0,        0]', ...
                                 obj.Lp(1)* [0.5,  sqrt(3)/2,   0]', ...
                        0,       obj.Lp(2)* [0,       0,        1]);
                case 'trigonal' %<111> is the 3fold axis
                    %    value =
                    c = cos(obj.Lp(1));
                    ty = sqrt((1. - c) / 6.);
                    tz = sqrt((1. + 2. * c) / 3.);
                    u = tz - 2. * sqrt(2.) * ty;
                    v = tz + sqrt(2.) * ty;
                    d = obj.Lp(1) / sqrt(3);
                    value = d*[u, v, v;
                               v, u, v;
                               v, v, u]';
                case 'tetragonal'
                    if strcmpi(obj.Centering, 'P')
                        value = cat2(obj.Lp(1)* [1, 0, 0], ...
                                     obj.Lp(1)* [0, 1, 0], ...
                                     obj.Lp(3)* [0, 0, 1]);
                    elseif strcmpi(obj.Centering, 'I')
                        a = obj.Lp(1) / 2.;
                        value = [a,  a, obj.Lp(3)/2.;
                                -a,  a, obj.Lp(3)/2.;
                                -a, -a, obj.Lp(3)/2.]';
                        value(:,:,2) = 0.5 * [1, -1, -1;
                            1,  1, -1;
                            1,  1,  1];
                    end
                case 'orthorhombic'
                    value = cat(2, obj.Lp(1)*[1 0 0]' , obj.Lp(2)*[0 1 0]' , obj.Lp(3)*[0 0 1]' );
                    % TODO
                    switch obj.Centering
                        case 'A' % 'B','C'  one-face centred, % TODO
                    %                     # base centered orthorhombic
                    %                     a = p[0]
                    %                     b = p[1]
                    %                     c = p[2]
                    %                     e1 = np.array([a / 2, b / 2, 0])
                    %                     e2 = np.array([-a / 2, b / 2, 0])
                    %                     e3 = np.array([0, 0, c])
                    %                     C = np.array([[0.5, -0.5, 0],
                    %                                   [0.5, 0.5, 0],
                    %                                   [0, 0, 1]])
                        case 'F'
                    %                     # face centered orthorhombic
                    %                     a = p[0]
                    %                     b = p[1]
                    %                     c = p[2]
                    %                     e1 = np.array([a / 2, b / 2, 0])
                    %                     e2 = np.array([0, b / 2, c / 2])
                    %                     e3 = np.array([a / 2, 0, c / 2])
                    %                     C = 0.5 * np.array([[1, 0, 1],
                    %                                         [1, 1, 0],
                    %                                         [0, 1, 1]])
                        case 'I'
                    %                     # body centered orthorhombic
                    %                     a = p[0]
                    %                     b = p[1]
                    %                     c = p[2]
                    %                     e1 = np.array([a / 2, b / 2, c / 2])
                    %                     e2 = np.array([-a / 2, b / 2, c / 2])
                    %                     e3 = np.array([-a / 2, -b / 2, c / 2])
                    %                     C = 0.5 * np.array([[1, -1, -1],
                    %                                         [1, 1, -1],
                    %                                         [1, 1, 1]])
                    end
% add other monoclinic settings
                 case 'monoclinic' % unique axis b
                     switch obj.Centering
                         case 'P'
                             a = obj.Lp(1);
                             b = obj.Lp(2);
                             c = obj.Lp(3);
                             beta = obj.Lp(3);
                             value = cat(2, [a, 0, 0]', [0, b, 0]', [c*cos(beta), 0, c*sin(beta)]');                             
                         case 'A' % 'B','C'  one-face centred, % TODO   # base centered monoclinic
                    %                     a = p[0]
                    %                     b = p[1]
                    %                     c = p[2]
                    %                     beta = radians(p[3])
                    %                     e1 = np.array([a / 2, b / 2, 0])
                    %                     e2 = np.array([-a / 2, b / 2, 0])
                    %                     e3 = np.array([c * np.cos(beta), 0, c * np.sin(beta)])
                    %                     C = np.array([[0.5, -0.5, 0],
                    %                                   [0.5, 0.5, 0],
                    %                                   [0, 0, 1]])
                     end
%                 case 'triclinic'
%                     a = obj.Lp(1);
%                     b = obj.Lp(2);
%                     c = obj.Lp(3);
%                     alpha = radians(p[3])
%                     beta = radians(p[4])
%                     gamma = radians(p[5])
%                     e1 = np.array([a, 0, 0])
%                     e2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
%                     e3 = np.array([c * np.cos(beta), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
%                         c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.cos(alpha) ** 2
%                     - np.cos(beta) ** 2 - np.cos(gamma) ** 2) / np.sin(gamma)])
                    %
                    %                 return np.array([e1, e2, e3]).T, la.inv(C) 
            end
            if ~exist('value', 'var')
                error('Bravais or centering type has not been set properly')
            end
            if size(value,3) == 1
                value(:,:,2) = eye(3);
            end
        end
        %------------------------------------------------------------------
        function base = get.C(obj)
            if size(obj.E,3) == 2
                base = obj.E(:,:,2)*obj.E(:,:,1);
            else
                base = obj.E(:,:,1);
            end
        end
        %------------------------------------------------------------------
        function metric = get.Cov_metric_C(obj)
            % given the basis vectors E of the primitive crystallographic lattice
            % calculate the corresponding covariant metric Z_i
            metric = calc_metric(obj.E);
        end
        %------------------------------------------------------------------
        function cont_metric = get.Contra_metric_C( obj )
            % given the basis vectors E of the primitive crystallographic lattice
            % calculate the corresponding Contravariant metric Z^i
            % mf = matrix_functions;
            cont_metric = inv(obj.Cov_metric_C);
        end
        %------------------------------------------------------------------
        function db = get.reciprocal_base(obj)
            % given the basis vectors E of the primitive crystallographic lattice
            % calculate its dual basis
            db = calc_dual_basis( obj.E );
        end
        %-------------------------------------------------------------------
        % The components of any vector referred, or the coordinates of any point
        % to the reciprocal basis represent the
        % Miller indices of a plane whose normal is along that vector, with the spacing
        % of the plane given by the inverse of the magnitude of that vector.
        %%%%%        
        % at the moment the primitive lattice is of no concern
        function prim = to_primitive_Miller( idx ) 
            % convert Miller Indices or a list of Miller indices stored in
            % an 3xN matrix from from conventional to primitive base
            % base to a primitive base
            for i = 1:size(idx,2)
                prim(:,i) = C * idx(:,i);
            end
        end
        %-------------------------------------------------------------------
        function miller = to_conventional_Miller( idx )
            % convert Miller indices or a list of Miller indices "idx"
            % stored in a 3xN matrix from primitive to conventional base
            % Miller indices are always integers!
            %mf = matrix_functions();
            for i = 1:size(idx,2)
                miller(:,i) = inv( obj.E(:,:,2) )* idx(:,i);
            end
        end
        %------------------------------------------------------------------
        function obj = set.density( obj, Mz) 
            % Mz... Matrix holding in the first row M_i
            % and in the second row z_i
            % M...atomic Mole Masses, z...number of atom species in the unit cell 
            volume = sqrt(det(Cov_metric));
            obj.density = (Mz(1,:).*Mz(2,:)) / volume * Avogadro; 
        end
        %------------------------------------------------------------------ 
        function vector = vec_from_coords(obj, coords )
            % assemble a lattice vector based on its components or a point
            % based on its coords
            vector = coords(1)*obj.C(:,1) + coords(2)*obj.C(:,2) + coords(3)*obj.C(:,3);
        end
        %------------------------------------------------------------------
        function length = vec_length( components )
            % given the components of a lattice vector calculates its length
            length = sqrt( components' * Cov_metric * components);
        end
        %------------------------------------------------------------------
        function distance = distance_atoms(coords1, coords2)
            X12 = coords1 - coords2;
            distance = vec_length( X12 );
        end
        %------------------------------------------------------------------
        function angle = bond_angle( coords1, coords2, coords3 )
            % given three atom positions. calculates the atomic angle so
            % thats it is measured at position of atom 1
            X12 = coords1 - coords2;
            X13 = coords1 - coords2;
            a = X12' * Cov_metric * X13; % Z_ij = U^i V^i
            angle = acos( a / (sqrt(a)* sqrt(X13' * Cov_metric * X12)));
        end
        %------------------------------------------------------------------
        function components = components_tensorcalculus( vec, basis)
            % given a vector and a co- or contravariant basis, this function yields the
            % co- or contravariant components of the vector. (co -> co, contra-> contra)
            % see Grinfeld - Intro to Tensorcalculus... p.61
            for i = 1:size(basis,2)
                components(i,1) = dot(vec, basis(:,i));
            end
        end
        %------------------------------------------------------------------
        % alternative solution of the above problem solving a coupled
        % system of equations...
        % function components = decomposition_wrt_to_basis( vec, basis )
        % % Given a vector (colon) and a basis (matrix with colon vectors)
        % % this function calculates the components of the vector w.r.t to the basis
        % % (for affine bases only, i.e. constant basis in space)
        % % See e.g. Grinfeld - Introduction to Tensor analysis calculus of moving surfaces - p61
        % components = zeros(3,1);
        % if isorthogonal(basis) % note: even simpler for an orthonormal basis
        %     for i=1:3
        %         components(i,1) = dot(vec , basis(:,i) ) / dot(basis(:,i),basis(:,i));
        %     end
        % else % for affine bases (i.e constant in space but non-orthogonal)
        %     vec_dot_e_i = zeros(3,1);
        %     for i=1:3
        %         vec_dot_e_i(i,1) = dot(vec , basis(:,i) );
        %     end
        %     metric = metric_from_basis( basis );
        %     components(i,1) = metric \ vec_dot_e_i;
        % end
        % end
    end % end methods
end % end class

