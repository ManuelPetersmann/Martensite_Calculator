classdef Composite_solution < IPS_solution
    % class for composite block solutions of two IPS_solutions
    % 
    % Note that the additive composition does only slightly change the IPS
    % condition but can highly influences the determinant!
    % class holds optimized fractions xi = [xi, 1-xi] of block solutions
    % for various cost/target functions
    % Fcompposite = xi * F1 + (1-xi)*F2
    
    properties (Access = public)
        %
        lath_id_pair; % then the lath solutions can be found
        % 
        % delta_hp_111
        %
        % rotation of average!
    end
    
    properties (Dependent) % all depend on lath_id_pair
        % The following cost-functions could be minimized w.r.t the phase fraction x, 1-x
        % and are stored in this class
        % 1 ) 'eps': the magnitude of the linearly mixed shape strain vector as a measure of
        % self accommodation.  % F_composite = I + eps * d_composite \otimes h_composite
        % min_norm2( eps1*d1 + eps2*d2)--> linear constrained optimization
        % 2 ) 'disg':  | F_composite -I |                   - displacement gradient
        % 3 ) 'gl':    | F_composite^T F_composite - I |    - c.f. Green-Lagrange (no rotation)
        %
        % -)F_composite should be mostly volumetric, i.e.  min|F_composite - % F^spherical| --> nonlinear constrained optimization
        % where the diagonal entries of F^spherical are determined from the volume
        % change (det(Bain)) - used e.g. in Qi2014 - though I think this is not reasonable
        
        x_eps; % optimized phase fractions of linear mixture
        eps_net; % minimum magnitude of shape vector by means of linear mixture of shape vectors of IPS solutions
        % F_comp_eps; % composite deformation gradient of above solution
        %
        x_dis; % minimum frobenius norm of displacement gradient for linear mixture Fcomp
        frob_opt_displacement_grad; % optimized phase fractions of linear mixture
        %
        x_gl; % % optimized phase fractions of linear mixture
        frob_opt_green_lagrange; % minimum frobenius norm of 2*green_lagrange tensor of linear mixture Fcomp        
    end 
    
    
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
%                 obj.eps_net = varargin{9};
%                 obj.x_eps   = varargin{10};
%                 %
%                 obj.x_dis = varargin{11};
%                 obj.frob_opt_displacement_grad = varargin{12};
%                 %
%                 obj.x_gl = varargin{13};
%                 obj.frob_opt_green_lagrange = varargin{14};

                obj.tolerances = containers.Map();
            end
        end
        %
        
        function x_eps = get.x_eps()
            % x = fmincon(@myfun,x0,A,b)
            [x_eps, eps_block] = mixture_vecs_lin_least_squares_opt( cat(2,sol1.eps_ips*sol1.d , sol2.eps_ips*sol2.d) );
            F_comp_eps = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_eps , 3,3) ;
            %x_eps
            %disp(['det(F_comp_eps) = ',num2str( det(F_comp_eps) ) ] );
            
            [x_dis, d_dis] = mixture_matrix_lin_least_squares_opt(cat(3,sol1.ST,sol2.ST));
            F_comp_dis = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_dis , 3,3) ;
            %x_dis
            %disp(['det(F_comp_dis) = ',num2str( det(F_comp_dis) ) ] );
            
            [x_gl,  d_gl] = frob_min_green_lagrange_composite_block( sol1.ST,sol2.ST, x_eps(1) );
            F_comp_gl = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_gl , 3,3) ;
            %x_gl
            %disp(['det(F_comp_gl) = ',num2str( det(F_comp_dis) ) ] );
            
            %         switch opt_func
            %             case 'eps'
            %                 F_comp = F_comp_eps;
            %             case 'disg'
            %                 F_comp = F_comp_dis;
            %             case 'gl'
            %                 F_comp = F_comp_gl;
            %         end
            
            % calculate composite IPS properties of all possibilites of optimized mixings
            try
                [y1_eps, y3_eps, d1_eps, d2_eps, h1_eps, h2_eps, Q1_eps, Q2_eps] = rank_one(F_comp_eps, I, lambda2_tol);
                % here a smaller tolerance value is used because the eigenvalue
                % changes slightly see introductory Note. If even this
                % tolerance is not reached the mixing is neglected
                
            catch %ME
                neglected_eigenvalue_change_during_mixing_opt_eps = neglected_eigenvalue_change_during_mixing_opt_eps + 2;
                %disp( ME );
                %[y1_eps, y3_eps, d1_eps, d2_eps, h1_eps, h2_eps, Q1_eps, Q2_eps] = deal(no_sol{:});
                %continue % do not safe solution if
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
                [y1_dis, y3_dis, d1_dis, d2_dis, h1_dis, h2_dis, Q1_dis, Q2_dis] = rank_one(F_comp_dis, I, lambda2_tol);
                % here a smaller tolerance value is used because the eigenvalue
                % changes slightly see introductory Note. If even this
                % tolerance is not reached the mixing is neglected
            catch %ME
                neglected_eigenvalue_change_during_mixing_opt_dis = neglected_eigenvalue_change_during_mixing_opt_dis + 2;
                %disp( ME );
                %[y1_dis, y3_dis, d1_dis, d2_dis, h1_dis, h2_dis, Q1_dis, Q2_dis] = deal(no_sol{:});
                %continue % do not safe solution if
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
                [y1_gl, y3_gl, d1_gl, d2_gl, h1_gl, h2_gl, Q1_gl, Q2_gl] = rank_one(F_comp_gl, I, lambda2_tol);
                % here a smaller tolerance value is used because the eigenvalue
                % changes slightly see introductory Note. If even this
                % tolerance is not reached the mixing is neglected
            catch %ME
                neglected_eigenvalue_change_during_mixing_opt_gl = neglected_eigenvalue_change_during_mixing_opt_gl + 2;
                %disp( ME );
                %[y1_gl, y3_gl, d1_gl, d2_gl, h1_gl, h2_gl, Q1_gl, Q2_gl] = deal(no_sol{:});
                %continue % do not safe solution if
            end
            
            
            no_sol = {0, 0, zeros(1,3), zeros(1,3), zeros(1,3), zeros(1,3), zeros(3), zeros(3) };
            
            y1 = [ y1_eps;  y1_dis;  y1_gl ], ...
                y3 = [ y3_eps;  y3_dis;  y3_gl ], ...
                d = [ d1_eps;  d2_eps;  d1_dis;  d2_dis; d1_gl;  d2_gl ], ...
                h = [ h1_eps;  h2_eps;  h1_dis;  h2_dis; h1_gl;  h2_gl ], cat(3,Q1_eps, Q2_eps, Q1_dis, Q2_dis, Q1_gl, Q2_gl)
                 
        end
     
    end % end methods
       
end 







