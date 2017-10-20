classdef Composite_solution < IPS_solution
    % class for composite block solutions of two IPS_solutions
    % 
    % Note that the additive composition does only slightly change the IPS
    % condition but can highly influence the determinant!
    % class holds optimized fractions xi = [xi, 1-xi] of block solutions
    % for various cost/target functions
    % Fcomp = xi * F1 + (1-xi)*F2
    
    properties (Access = public)
        %
        lath_id_pair;
        %
        eps_net; % minimum magnitude of shape vector by means of linear mixture of shape vectors of IPS solutions
        x_eps; % optimized phase fractions of linear mixture
        F_comp_eps; % composite deformation gradient of above solution
        %
        x_dis; % minimum frobenius norm of displacement gradient for linear mixture Fcomp
        frob_opt_displacement_grad; % optimized phase fractions of linear mixture
        %
        x_gl; % % optimized phase fractions of linear mixture
        frob_opt_green_lagrange; % minimum frobenius norm of 2*green_lagrange tensor of linear mixture Fcomp        
    end % end of properties
    
    methods
        % constructor
        function obj = Composite_solution( varargin ) % IPS_Solution{ F, G, id, eps_0, d, h, Q, LT } + Composite solution variables
            %
            if nargin == 0 % =  if isempty(varargin), return;end
                % no argument constructor
                super_args = {};
            end
            if nargin > 8
                super_args = { varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8} };
            end
            %
            obj = obj@IPS_solution( super_args{:} ); % actually only needs: F, G, id, eps_0, d, h, Q, LT
            %
            if nargin > 8
                obj.eps_net = varargin{1,9};
                obj.x_eps   = varargin{1,10};
                obj.F_comp_eps = varargin{1,11};
                %
                obj.x_dis = varargin{1,12};
                obj.frob_opt_displacement_grad = varargin{1,13};
                %
                                obj.x_gl = varargin{1,15};
                obj.frob_opt_green_lagrange = varargin{1,14};

            end
        end
    end % end methods
       
end 







