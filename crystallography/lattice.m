classdef Lattice  
    % This class is for the specification of lattice type --> symmetry
    % and lattice parameters
    
    properties
        Bravais_type % a string specifying the Bravais lattice type e.g. 'cubic'
        Centering % a string specifying the exact type e.g. simple 'P', face-centered 'I' etc.
        % Lattice_parameters
        Lp = {1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
            'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
        % cell indexing with ()
        % content indexing with {}
        % acces with Lattice_params(1:2) -> 1, 'a'
        C % must be user specified because of different possibilities
        % conventional basis c_i stored as C = [c1;c2;c3], (such that as
        % many angles as possible amout to 90°. e.g. 
        % primitive lattice --> conventinal lattice = primitive (e.g. primitive cubic),
        % otherwise "centered" e.g. face centered cubic.
    end
    properties (Dependent, Access = public)
        Point_group  % 3x3xN matrix containing the point group symmetry matrices
        E
        % primitive (crystallographic) basis vectors e_i stored as E =
        % [e1;e2;e3], row vectors
        % i.e. lattice vectors with integer coefficients t_i: {t_vec = t1*a1 + t2*a2 + t3*a3}.
        % this means that the a_i are translation vectors!
        Metric
    end
    
    methods % in matlab no special or hidden class object is passed to all methods (only explicitly)
        % obj.somemethod(arg) == somemethod(obj, arg)
        % .set functions control input 
        
        % a = obj.('Propertyname') --> use this to pass property names as
        % arguments
        
        % constructor
        function obj = Lattice()
        end
        %-------------------------------------------------------------------
        % outside-class function-signatures:
        pointgroup = pointgroups( Bravais_type )
        
        %-------------------------------------------------------------------
        %TODO write function or class for user input e.g.
        %bravais_type = input('choose bravais lattice') %later maybe make a
        %GUI for that#
        %-------------------------------------------------------------------
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
        %-------------------------------------------------------------------
        function obj = set.Centering(obj, input )
            Centering_strings = ['P'; 'C'; 'I'; 'F'];
            if ismember( input, Centering_strings)
                obj.Centering = input;
            else
                error('Invalid Centering')
            end
        end
        %-------------------------------------------------------------------
%         function obj = set.Lp(obj)
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
%             end                         
%         end    
        %-------------------------------------------------------------------
        %function obj = set.C(obj)
        %     
        %end
        %-------------------------------------------------------------------
        
        function point_group = get.Point_group( obj ) 
            pg = pointgroups;
            point_group = pg.determine_pointgroup( obj.Bravais_type );
        end
        %-------------------------------------------------------------------
        % TODO
        %         function obj = get.Metric(obj)
        %         %TODO given the basis vectors calculate the metric
        %         end
        %-------------------------------------------------------------------
        %         function [bool] = is_isometry( M )
        %             Given an affine Mapping M, calculates wheter for the given
        %             lattice the mapping is an isometry
        %             TODO calc metric
        %             if G == W'*G*W --> isometry
        %         end
        %-------------------------------------------------------------------
        function obj = get.E(obj)
        % generate the base matrix for primitive cell of 3D Bravais lattices
            switch obj.Bravais_type
                case 'cubic P'
                case 'cubic'
                    switch obj.Centering
                        case 'P' % simple
                            obj.E = obj.Lp{1} * eye(3);
                        case 'F' % Face centered fcc
                            obj.E = 0.5 * obj.Lp{1} * eye(3);
                        case 'I' % Body centered bcc
                            obj.E = 0.5 * obj.Lp{1} * [1, -1, -1;
                                                       1,  1, -1;
                                                       1,  1,  1];
                    end
                case 'hexagonal'
                    obj.E = [obj.Lp{1}, obj.Lp{1}*0.5           0;
                                0,      obj.Lp{1}*sqrt(3)/2     0;
                                0,          0,          obj.Lp{2}];
                %case 'trigonal' %<111> is the 3fold axis
                %    obj.E = 
%                     c = np.cos(p[1] * np.pi / 180)
%                     a = p[0]
%                     ty = np.sqrt((1 - c) / 6)
%                     tz = np.sqrt((1 + 2 * c) / 3)
%                     u = tz - 2 * np.sqrt(2) * ty
%                     v = tz + np.sqrt(2) * ty
%                     e1 = a / np.sqrt(3) * np.array([u, v, v])
%                     e2 = a / np.sqrt(3) * np.array([v, u, v])
%                     e3 = a / np.sqrt(3) * np.array([v, v, u])
                %case
            end
            
        end
        %-------------------------------------------------------------------          
    end % end of methods    
end % end of class


%                 elif n == 5:
%                     # trigonal
%                     

%                 elif n == 6:
%                     # simple tetragonal
%                     a = p[0]
%                     c = p[1]
%                     e1 = a * np.array([1, 0, 0])
%                     e2 = a * np.array([0, 1, 0])
%                     e3 = c * np.array([0, 0, 1])
%                 elif n == 7:
%                     # body centered tetragonal
%                     a = p[0]
%                     c = p[1]
%                     e1 = (a / 2) * np.array([1, 1, c / a])
%                     e2 = (a / 2) * np.array([-1, 1, c / a])
%                     e3 = (a / 2) * np.array([-1, -1, c / a])
%                     C = 0.5 * np.array([[1, -1, -1],
%                                         [1, 1, -1],
%                                         [1, 1, 1]])
%                 elif n == 8:
%                     # simple orthorhombic
%                     a = p[0]
%                     b = p[1]
%                     c = p[2]
%                     e1 = np.array([a, 0, 0])
%                     e2 = np.array([0, b, 0])
%                     e3 = np.array([0, 0, c])
%                 elif n == 9:
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
%                 elif n == 10:
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
%                 elif n == 11:
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
%                 elif n == 12:
%                     # monoclinic unique axis b
%                     a = p[0]
%                     b = p[1]
%                     c = p[2]
%                     beta = radians(p[3])
%                     e1 = np.array([a, 0, 0])
%                     e2 = np.array([0, b, 0])
%                     e3 = np.array([c * np.cos(beta), 0, c * np.sin(beta)])
%                 elif n == 13:
%                     # base centered monoclinic
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
%                 elif n == 14:
%                     # triclinic
%                     a = p[0]
%                     b = p[1]
%                     c = p[2]
%                     alpha = radians(p[3])
%                     beta = radians(p[4])
%                     gamma = radians(p[5])
%                     e1 = np.array([a, 0, 0])
%                     e2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
%                     e3 = np.array([c * np.cos(beta), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
%                                    c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.cos(alpha) ** 2
%                                                - np.cos(beta) ** 2 - np.cos(gamma) ** 2) / np.sin(gamma)])
%         
%                 return np.array([e1, e2, e3]).T, la.inv(C)

        
        
        
