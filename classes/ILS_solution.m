classdef ILS_solution 
    % ILS_solution
    
    properties (Access = public)
        k; % invariant line vector
        ST; % (ST= R*B*S2*S1) u = u   % shape transformation for invariant line
        LT; % lattice transformation - RB
        R_Bain; % rotation necessary to bring  B*u - back to u (same direction but in general stretched!)
        %
        added_props; % = containers.Map();
        id; % to know after sorting which solutions came in pairs initially 1-2, 3-4, 5-6... and for block solutions
        %
        slip; % Slip_systems()
        shear_increments;
    end
    properties (Dependent)
        lambda2_IPS_to_one; % % characterises how close the ILS is to an IPS
        R_inclusion;
        axis_angle_rotvec_inclusion; % returns 1x4 vector of rotation axis [1:3] and angle in degree [4]
        rotangle_inclusion; % same as above but specifically used for reduction of solutions
        % for this term a cosserat like contribution to the strain energy could be written,
        % however preferable it should be vanishingly small in reality!
        % here the rotation due to the shear is incorporated!!!
        R_lattice; % rotation matrix (for invariant line), without rotation due to shear! --> lattice rotation
        frob_green_lagrange;
        frob_displacement_grad;
    end
    
    methods
        % constructor: The constructor can return only a single argument.
        function obj = ILS_solution( varargin ) % F1, F2, id, eps, d, h, Q, LT)
            if isempty(varargin)
                %return; % no argument constructor
            else
            %obj.added_props = containers.Map; % assignment is done in class Solution_array
            obj.u = varargin{1};
            obj.ST = varargin{2};
            obj.LT = varargin{3}; % the lattice transformation = Q*Bain is given here directly 
            % so that the class does not need the Bain strain as property (only needed for this)
            obj.R_Bain = varargin{4};
            end
        end
        %
        %% get functions
        function R = get.R_lattice( obj )
            [~,R] = polardecomposition( obj.LT );
        end
        %
        function R = get.R_inclusion( obj )
            [~,R] = polardecomposition( obj.ST );
        end
        %
        function l2 = get.lambda2_IPS_to_one( obj )
            [ ~, y2] = sorted_eig_vals_and_vecs( obj.ST' * obj.ST );
            l2 = abs(y2 -1.);
        end
        %
        function frobgl = get.frob_green_lagrange(obj)
            gl = obj.ST' * obj.ST - eye(3);
            frobgl = trace(gl'*gl);        
        end
        %
        function frobdis = get.frob_displacement_grad(obj)
            dis = obj.ST - eye(3);
            frobdis = trace(dis' * dis);
        end
        %
        function vec4 = get.axis_angle_rotvec_inclusion( obj )
            [~,Q] = polardecomposition( obj.ST );
            [angle, axis] = rotmat_to_axis_angle( Q ); %vrrotmat2vec( Q );
            vec4 = cat(1,axis,angle);
        end
        %
        function abs_angle = get.rotangle_inclusion( obj )
            [~,Q] = polardecomposition( obj.ST );
            abs_angle = acosd( (trace(Q)-1.) / 2.);
        end

    end % methdos


    
end % class



