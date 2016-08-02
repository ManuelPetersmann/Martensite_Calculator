classdef IPS_solution < dynamicprops
    % SOLUTION - Baseclass for solutions of the IPS equation: Q*F - G = eps0 * (a \dyad n)
    
    properties (Access = public)
        F = zeros(3); % e.g. modified Bain strain BS or 0.5*( R*S + inverse(R)*S_mirror ) * B as in Qi2013
        G = zeros(3); % this should be the reference deformation, e.g. Identiy for austenite, or fixed U_i for twins
        id = 0;
        eps = 0.; % strain: lambda_3 - lambda_1 or sqrt respectively depending on convention
        n = [0. 0. 0.]'; % normal vector to habit plane (unit vector)
        a = [0. 0. 0.]'; % shear direction of transformation (unit vector)
        Q = zeros(3); % rotation matrix (for invariant planar match between domains of homogeneous deformation F and G
        LT = zeros(3); % calculation of Lattice-Transformation (A_L in Qi2014 Paper)
        %LT; % = RB = ...R von IPS condition
    end
    properties (Dependent)
        % shape transformation (A_D in Qi2014 Paper)
        ST;  % on side of homogeneous deformation F:  QF - G = eps_0* ( a \dyad n )
    end
    
    methods
        % constructor: The constructor can return only a single argument.
        function obj = IPS_solution( varargin ) % F, G, id, eps, a, n, Q, LT)
            if isempty(varargin)
                return; % no argument constructor
            end
            %
            if nargin == 8 % constructor from solution arguments
                obj.F = varargin{1};
                obj.G  = varargin{2};
                obj.id = varargin{3};
                obj.eps= varargin{4};
                obj.n  = varargin{5};
                obj.a  = varargin{6};
                obj.Q  = varargin{7};
                obj.LT = varargin{8}; % the lattice transformation = Q*Bain is given here directly so that the class does not need the Bain strain as property (only needed for this)
            end
        end
        
        %%
        function shape_transformation = get.ST( obj )
            shape_transformation = obj.Q * obj.F; % = G + eps_0 ( a \dyad n )
        end
        %
        %         function lattice_transformation = get.LT( obj )
        %             lattice_transformation = obj.Q * obj.U; % = G + eps_0 ( a \dyad n )      %%%%%%%%% Problem like this is that the bain strain would occur in two classes...
        %         end                                                                              %%%%%%%%% Or, if this class is derived from martensite - the object array contains to much data...
        
        
    end % methdos
end % class



