classdef IPS_solution 
    % IPS_solution - Baseclass for solutions of the IPS equation: F1 - F2 = eps0 * (d \dyad h)
    
    properties (Access = public)
        F1; % F1=B*S2*S1 for double shear or 0.5*( R*S + inverse(R)*S_mirror ) * B as in Qi2013
        F2; % F2=U2 R2 this should be the reference deformation, e.g. Identiy for austenite, or fixed U_i for twins
        %
        % PET 10.10.17: replaced 'id' and 'eps_ips' wit lambda_1 and lambda_2,
        % made 'eps_ips' = sqrt(lambda_3) - sqrt(lambda_1) (or sqrt respectively
        % depending on convention) dependent property 
        lambda_1; % smallest eigenvalue (deformation)
        %lambda_2 = 1 !-> IPS
        lambda_3; % largest  eigenvalue (deformation)
        %
        h; % normal vector to habit plane (unit vector)
        d; % shear direction of transformation (unit vector)
        Q; % rotation matrix (for invariant planar match between domains of homogeneous deformation F and G
        LT; % calculation of Lattice-Transformation (A_L in Qi2014 Paper)
        %LT; % = RB = ...R von IPS condition
        %
        added_props; % PET: 19.10.17 | = containers.Map();
        id; % to know after sorting which solutions came in pairs initially 1-2, 3-4, 5-6... and for block solutions
        % could then be incremented in the constructor...
    end
    properties (Dependent)
        % shape transformation (A_D in Qi2014 Paper)
        eps_ips;
        ST;  % on side of homogeneous deformation F:  QF - G = eps_0* ( a \dyad n )
        % note those are the eigenvectors of F1'*F1 hence they do not
        % include any rotation!!! the two lines below have therefore been commented - not valid!
        %e1; %dir_of_largest_def; %  = eigenvector corresponding to lambda_3 - should be close to the habit plane vector (|| shortest dimension 'c' of lath)
        %e3; %dir_of_smallest_def; % = eigenvector corresponding to lambda_1 - should be close to largest dimension 'a' of lath
        frob_green_lagrange;
        frob_displacement_grad;
        axis_angle_rotvec_inclusion; % returns 1x4 vector of rotation axis [1:3] and angle in degree [4]
        % for this term a cosserat like contribution to the strain energy could be written,
        % however preferable it should be vanishingly small in reality!
        %
        % angle_e1_to_cpdir; % angle between smallest deformation and close packed_direction; %to_Invariant line strain _KS_NW_dir_aust;
        % added via the class Solution_array_dynamically (needs information
        % from austenite object...
    end
    
    methods
        % constructor: The constructor can return only a single argument.
        function obj = IPS_solution( varargin ) % F1, F2, id, eps, d, h, Q, LT)
            if isempty(varargin)
                return; % no argument constructor
            end
            %
            if nargin == 8 % constructor from solution arguments
                %
                % obj.added_props = containers.Map; % assignment is done in
                % class Solution_array
                % obj.id = count_instance_id( 1 ); % add plus one to id for newly created IPS_solution
                % this is better done directly in the middle_eigenvalue_mod function
                obj.F1 = varargin{1};
                obj.F2  = varargin{2};
                obj.lambda_1 = varargin{3};
                obj.lambda_3 = varargin{4};
                obj.d  = varargin{5};
                obj.h  = varargin{6};
                obj.Q  = varargin{7};
                obj.LT = varargin{8}; % the lattice transformation = Q*Bain is given here directly so that the class does not need the
                                      % Bain strain as property (only needed for this)
            end
        end
        %
        %% get functions
        function shape_transformation = get.ST( obj )
            shape_transformation = obj.Q * obj.F1; % = G + eps_0 ( a \dyad n )
        end
        %
        function eps_ips = get.eps_ips( obj )
            eps_ips = sqrt(obj.lambda_3) - sqrt(obj.lambda_1);
        end
        %
        function frobgl = get.frob_green_lagrange(obj)
            frobgl = trace(obj.ST' * obj.ST - eye(3));        
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
%         function smallest_eigenvector = get.e1( obj )
%             [~,~,~,smallest_eigenvector] = sorted_eig_vals_and_vecs( obj.F1'*obj.F1 ); % [ y1, y2, y3, e1, e2, e3] = 
%         end
%         %
%         function largest_eigenvector = get.e3( obj )
%             [~,~,~,~,~,largest_eigenvector] = sorted_eig_vals_and_vecs( obj.F1'*obj.F1 ); % [ y1, y2, y3, e1, e2, e3] =
%         end
%%  is not done like this, but added via the class Solution_array_dynamically   
%        function ang = get.angle_e1_to_cpdir(obj)
%            %if isprop(obj,'closest_KS') % generally cp-direction
%            min_misorientation( cpps_gamma, obj.e1 );
%            %end
%        end
        %         function lattice_transformation = get.LT( obj )
        %             lattice_transformation = obj.Q * obj.U; % = G + eps_0 ( a \dyad n )   %%%%%%%%% Problem like this is that the bain strain would occur in two classes...
        %         end                                                     %%%%%%%%% Or, if this class is derived from martensite - the object array contains to much data...       
        
    end % methdos

% If i do it like this - how / when do i set it back to zero?
%     methods (Static, Access = private)
%         function new_id = count_instance_id(increment)
%             persistent VALUE
%             if isempty(VALUE)
%                 VALUE = 0;
%             end
%             if nargin > 0
%                 VALUE = VALUE + increment;
%             end
%             new_id = VALUE;
%         end
%     end
    
end % class



