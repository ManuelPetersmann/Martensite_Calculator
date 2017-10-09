classdef IPS_solution < dynamicprops
    % needs to be derived from dynamicprops because slip_solution can dynamically add e.g. OR directions!
    % IPS_solution - Baseclass for solutions of the IPS equation: F1 - F2 = eps0 * (d \dyad h)
    
    properties (Access = public)
        F1; % F1=U1 R1 e.g. modified Bain strain BS or 0.5*( R*S + inverse(R)*S_mirror ) * B as in Qi2013
        F2; % F2=U2 R2 this should be the reference deformation, e.g. Identiy for austenite, or fixed U_i for twins
        id;
        eps_ips; % strain: lambda_3 - lambda_1 or sqrt respectively depending on convention
        h; % normal vector to habit plane (unit vector)
        d; % shear direction of transformation (unit vector)
        Q; % rotation matrix (for invariant planar match between domains of homogeneous deformation F and G
        LT; % calculation of Lattice-Transformation (A_L in Qi2014 Paper)
        %LT; % = RB = ...R von IPS condition
    end
    properties (Dependent)
        % shape transformation (A_D in Qi2014 Paper)
        ST;  % on side of homogeneous deformation F:  QF - G = eps_0* ( a \dyad n )
        dir_of_largest_def; %  = eigenvector corresponding to the maximum eigenvalue - should be close to the habit plane vector
        dir_of_smallest_def; % = eigenvector corresponding to the minimum eigenvalue - should be close to long direction (a) of lath
        frob_green_lagrange;
        frob_displacement_grad;
        angle_smallest_def_to_close_packed_direction; %to_ILS_KS_NW_dir_aust;
        axis_angle_rotvec_inclusion; % returns 1x4 vector of rotation axis [1:3] and angle in degree [4]
        % for this term a cosserat like contribution to the strain energy could be written,
        % however preferable it should be vanishingly small in reality!
    end
    
    methods
        % constructor: The constructor can return only a single argument.
        function obj = IPS_solution( varargin ) % F1, F2, id, eps, d, h, Q, LT)
            if isempty(varargin)
                return; % no argument constructor
            end
            %
            if nargin == 8 % constructor from solution arguments
                obj.F1 = varargin{1};
                obj.F2  = varargin{2};
                obj.id = varargin{3};
                obj.eps_ips= varargin{4};
                obj.d  = varargin{5};
                obj.h  = varargin{6};
                obj.Q  = varargin{7};
                obj.LT = varargin{8}; % the lattice transformation = Q*Bain is given here directly so that the class does not need the Bain strain as property (only needed for this)
            end
        end
        %% get functions
        function shape_transformation = get.ST( obj )
            shape_transformation = obj.Q * obj.F1; % = G + eps_0 ( a \dyad n )
        end
        %
        function smallest_eigenvector = get.dir_of_smallest_def( obj )
            [~,~,~,smallest_eigenvector] = sorted_eig_vals_and_vecs( obj.F1 ); % [ y1, y2, y3, e1, e2, e3] = 
        end
        %
        function largest_eigenvector = get.dir_of_largest_def( obj )
            [~,~,~,~,~,largest_eigenvector] = sorted_eig_vals_and_vecs( obj.F1 ); % [ y1, y2, y3, e1, e2, e3] =
        end
        %
        function frobgl = get.frob_green_lagrange(obj)
            frobgl = trace(obj.ST^T * obj.ST - eye(3));        
        end
        %
        function frobdis = get.frob_displacement_grad(obj)
            frobdis = trace(obj.ST - eye(3));
        end
        %
        function vec4 = get.axis_angle_rotvec_inclusion( obj )
            [~,R] = polardecomposition( obj.ST );
            vec4 = vrrotmat2vec( R );
            vec4(4) = rad2deg( vec4(4) );
        end
        %
        function ang = get.angle_smallest_def_to_close_packed_direction(obj)
            %if isprop(obj,'closest_KS')
                ang = get_angle( obj.dir_of_smallest_def,obj.nearest_cp_direction );
            %end
        end
        %         function lattice_transformation = get.LT( obj )
        %             lattice_transformation = obj.Q * obj.U; % = G + eps_0 ( a \dyad n )   %%%%%%%%% Problem like this is that the bain strain would occur in two classes...
        %         end                                                                       %%%%%%%%% Or, if this class is derived from martensite - the object array contains to much data...       
        
    end % methdos
end % class



