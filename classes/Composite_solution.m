classdef Composite_solution < IPS_solution
    % class for composite block solutions of two IPS_solutions
    % 
    % Note that the additive composition does only slightly change the IPS
    % condition but can highly influences the determinant!
    % class holds optimized fractions xi = [xi, 1-xi] of block solutions
    % for various cost/target functions
    % Fcompposite = xi * F1 + (1-xi)*F2
    
    properties (Access = public)
        s; %lath_solution_pair; % type (e.g. IPS_solution or Slip_solution) is set in the constructor
    end
    
    properties (Dependent) % all depend on lath_id_pair
             
        % The following cost-functions could be minimized w.r.t the phase fraction x, 1-x
        % and are stored in this class:
        %
        % 1 ) 'eps': the magnitude of the linearly mixed shape strain vector as a measure of
        % self accommodation.  % F_composite = I + eps * d_composite \otimes h_composite
        % min_norm2( eps1*d1 + eps2*d2)--> linear constrained optimization
        shape_vec_opt; % =[ x_eps, eps_net]
        % x_eps; % [x, 1-x] optimized phase fractions of linear mixture
        % eps_net; % minimum magnitude of shape vector by means of linear mixture of shape vectors of IPS solutions
        
        % 2 ) 'disg':  | F_composite -I |                   - displacement gradient
        disg_opt; % = [x_dis, frob_opt_displacement_grad]
        % x_dis; % minimum frobenius norm of displacement gradient for linear mixture Fcomp
        % frob_opt_displacement_grad; % optimized phase fractions of linear mixture
        
        % 3 ) 'gl':    | F_composite^T F_composite - I |    - c.f. Green-Lagrange (no rotation)
        gl_opt; % = [x_gl, frob_opt_green_lagrange]
        % x_gl; % % optimized phase fractions of linear mixture
        % frob_opt_green_lagrange; % minimum frobenius norm of 2*green_lagrange tensor of linear mixture Fcomp  
        
        % -)F_composite should be mostly volumetric, i.e.  min|F_composite - % F^spherical| --> nonlinear constrained optimization
        % where the diagonal entries of F^spherical are determined from the volume
        % change (det(Bain)) - used e.g. in Qi2014 - though I think this is not reasonable

        
        %% Alternatively any phase fraction between [0,1] can be given and the average HP parameters are calculated
           delta_hp;  % get function takes cryst family, e.g. {111} and gives min misorientation to it.
           %rotation of average --> in IPS class!
    end 
    
    
    methods
        %% no constructor needed apart from a default one without arguments for now
        % everything is assigned in the set.lath_solution_pair function
        
        function obj = set.lath_solution_pair(obj, lath_sols)
            lath_class_name = class( lath_sols(1) ); % output string
            switch lath_class_name
                case 'IPS_Solution'
                    obj.lath_solution_pair = IPS_Solution();
                case 'Slip_solution'
                    obj.lath_solution_pair = Slip_solution();
            end
            obj.lath_solution_pair = lath_sols; 
        end
          
%         function obj = set.lath_id_pair(obj, ids)
%             if length(ids) ~= 2
%                 error('excatly two ids are needed');
%             end
%             % integer check
%             if ~all(mod(ids,1)== 0)
%                 error('lath ids must be integers!');
%             end
%         end
                          
%         % function to get full lath solution from lath id given the array lath_sols
%                 function laths = get.lath_solution_pair( lath_sols )
%                     % lath_sols must be an input!                    
%                     indices = find(  ( lath_sols.array.id == obj.lath_id_pair(1) | lath_sols.array.id == obj.lath_id_pair(2) )  );
%                     laths(1) = lath_sols.array( indices(1) );
%                     laths(2) = lath_sols.array( indices(2) );
%                 end
            
            
            
%             obj = obj@IPS_solution( super_args{:} ); % actually only needs: F, G, id, eps_0, d, h, Q, LT
%             %
%             if nargin > 8
%                                 obj.eps_net = varargin{9};
%                                 obj.x_eps   = varargin{10};
%                                 %
%                                 obj.x_dis = varargin{11};
%                                 obj.frob_opt_displacement_grad = varargin{12};
%                                 %
%                                 obj.x_gl = varargin{13};
%                                 obj.frob_opt_green_lagrange = varargin{14};
%                 
%                 obj.tolerances = containers.Map();
%             end
%         end


                % IPS_Solution{ F, G, id, eps_0, d, h, Q, LT } + Composite solution variables
                
                %                 super_args = {F1, F2, Bain, tolerance};
                %                 { varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8} };
                
                %                 block_solutions.array( isol )   =  Composite_solution(F1, I, y1, y3, d, h, Q, ...
                %                     zeros(3), eps_block, x_eps, F1, x_dis, d_dis, x_gl, d_gl );
                

        

        %% deviation of average block habit plane from given crystallographic set/family
        function dhp = get.delta_hp( x, set )
            set
            [y1, y3, d1, d2, h1, h2, Q1, Q2] = rank_one(Fc, I, lambda2_tol, false); % last 'false' is that no lambda_2_warning occurs
            min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h2) > block_hp_cp_aust_tol
        end
        
        
        % x = fmincon(@myfun,x0,A,b)
        function svo = get.shape_vec_opt()
            s1 = lath_solution_pair(1);
            s2 = lath_solution_pair(2);
            [x_eps, eps_block] = mixture_vecs_lin_least_squares_opt( cat(2,s1.eps_ips*s1.d , s2.eps_ips*s2.d) );
            %F_comp_eps = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_eps , 3,3) ;
            %x_eps
            %disp(['det(F_comp_eps) = ',num2str( det(F_comp_eps) ) ] );
        end
            
                function dgo = get.disg_opt()
            [x_dis, d_dis] = mixture_matrix_lin_least_squares_opt(cat(3,sol1.ST,sol2.ST));
            %F_comp_dis = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_dis , 3,3) ;
            %x_dis
            %disp(['det(F_comp_dis) = ',num2str( det(F_comp_dis) ) ] );
                end
            
                    function glo = get.gl_opt()
            [x_gl,  d_gl] = frob_min_green_lagrange_composite_block( sol1.ST,sol2.ST, x_eps(1) );
            %F_comp_gl = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_gl , 3,3) ;
            %x_gl
            %disp(['det(F_comp_gl) = ',num2str( det(F_comp_dis) ) ] );
                    end
                    
            
            %         switch opt_func
            %             case 'eps'
            %                 F_comp = F_comp_eps;
            %             case 'disg'
            %                 F_comp = F_comp_dis;
            %             case 'gl'
            %                 F_comp = F_comp_gl;
            %         end
            
            % calculate composite IPS properties of all possibilites of optimized mixings
            
            [y1_eps, y3_eps, d1_eps, d2_eps, h1_eps, h2_eps, Q1_eps, Q2_eps] = rank_one(F_comp_eps, I, lambda2_tol);
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [y1_dis, y3_dis, d1_dis, d2_dis, h1_dis, h2_dis, Q1_dis, Q2_dis] = rank_one(F_comp_dis, I, lambda2_tol);
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [y1_gl, y3_gl, d1_gl, d2_gl, h1_gl, h2_gl, Q1_gl, Q2_gl] = rank_one(F_comp_gl, I, lambda2_tol);
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
            
            
            
            %  no_sol = {0, 0, zeros(1,3), zeros(1,3), zeros(1,3), zeros(1,3), zeros(3), zeros(3) };
            %
            % y1 = [ y1_eps;  y1_dis;  y1_gl ], ...
            % y3 = [ y3_eps;  y3_dis;  y3_gl ], ...
            % d = [ d1_eps;  d2_eps;  d1_dis;  d2_dis; d1_gl;  d2_gl ], ...
            % h = [ h1_eps;  h2_eps;  h1_dis;  h2_dis; h1_gl;  h2_gl ], cat(3,Q1_eps, Q2_eps, Q1_dis, Q2_dis, Q1_gl, Q2_gl)
            
        end
        
    end % end methods
       
end 







