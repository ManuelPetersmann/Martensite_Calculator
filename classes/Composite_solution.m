classdef Composite_solution %< IPS_solution
    % class for composite block solutions of two IPS_solutions
    %
    % Note that the additive composition does only slightly change the IPS
    % condition but can highly influences the determinant!
    % class holds optimized fractions xi = [xi, 1-xi] of block solutions
    % for various cost/target functions
    % Fcompposite = xi * F1 + (1-xi)*F2
    
    properties %(Access = public)
        lath_solution_pair; % type (e.g. IPS_solution or Slip_solution) is set in the constructor
        Fc05; % linear mixture of deformations at xi=0.5
    end
    
    properties (Dependent) % all depend on lath_id_pair
        
        % The following cost-functions could be minimized w.r.t the phase fraction x, 1-x
        % and are stored in this class:
        
        % 1 ) 'eps': the magnitude of the linearly mixed shape strain vector as a measure of
        % self accommodation.  % F_composite = I + eps * d_composite \otimes h_composite
        % min_norm2( eps1*d1 + eps2*d2)--> linear constrained optimization
        shape_vec_opt; % =[ x_eps, eps_net, delta_hps_svo, rot_svo]
        % x_eps - [x, 1-x] optimized phase fractions of linear mixture
        % eps_net -  minimum magnitude of shape vector by means of linear mixture of shape vectors of IPS solutions
        % svo_delta_hps  -  get function takes cryst family, e.g. {111} and gives min misorientation to it.
        % svo_rot -  rotation of optimized average 
        
        
        % 2 ) 'disg':  | F_composite -I |                   - displacement gradient
        disg_opt; % = [x_dis, frob_opt_displacement_grad, disg_delta_hps, rot_svo]
        % x_dis  -   minimum frobenius norm of displacement gradient for linear mixture Fcomp
        % frob_opt_displacement_grad  -  optimized phase fractions of linear mixture
        % disg_delta_hps  -  get function takes cryst family, e.g. {111} and gives min misorientation to it.
        % disg_rot -  rotation of optimized average 
        
        
        % 3 ) 'gl':    | F_composite^T F_composite - I |    - c.f. Green-Lagrange (no rotation)
        gl_opt; % = [x_gl, frob_opt_green_lagrange, gl_delta_hps, rot_svo]
        % x_gl - optimized phase fractions of linear mixture
        % frob_opt_green_lagrange -  minimum frobenius norm of 2*green_lagrange tensor of linear mixture Fcomp
        % gl_delta_hps - get function takes cryst family, e.g. {111} and gives min misorientation to it.
        % gl_rot - rotation of optimized average 
        
        
        % -)F_composite should be mostly volumetric, i.e.  min|F_composite - % F^spherical| --> nonlinear constrained optimization
        % where the diagonal entries of F^spherical are determined from the volume
        % change (det(Bain)) - used e.g. in Qi2014 - though I think this is not reasonable
    end
    
    
    methods
        
        function obj = set.lath_solution_pair(obj, sols) % Composite_solution( sol1, sol2) 
            lath_class_name = class( sols(1) ); % output string
            switch lath_class_name
                case 'IPS_Solution'
                    obj.lath_solution_pair = IPS_Solution();
                case 'Slip_solution'
                    obj.lath_solution_pair = Slip_solution();
            end
            obj.lath_solution_pair = sols;
            % NECESSARY ??? - eps_ips almost doesn't change...
            % obj = obj@IPS_solution( super_args{lath_sols, Bain, tolerance} ); 
            % obj.tolerances = containers.Map();
        end  
        
%         function Fc05 = get.Fc05( obj )
%          Fc05 = linmix2( 0.5, obj.lath_solution_pair(1).ST, obj.lath_solution_pair(2).ST );
%         end
        %% phase fraction (x) optimization methods
        % x = fmincon(@myfun,x0,A,b)
        function svo = get.shape_vec_opt( obj )
            [x_eps, eps_block] = mixture_vecs_lin_least_squares_opt( ...
                cat(2,obj.lath_solution_pair(1).eps_ips*obj.lath_solution_pair(1).d , obj.lath_solution_pair(2).eps_ips*obj.lath_solution_pair(2).d) );
            %
            Fc = linmix2( x_eps, obj.lath_solution_pair(1).ST, obj.lath_solution_pair(2).ST );
            svo_delta_hps = delta_hps( Fc );
            svo_rot = axis_angle_rotvec_inclusion( Fc );
            svo = [ x_eps', eps_block, svo_delta_hps, svo_rot];
        end
        
        function dgo = get.disg_opt( obj )
            [x_dis, d_dis] = mixture_matrix_lin_least_squares_opt( ...
                cat(3,obj.lath_solution_pair(1).ST, obj.lath_solution_pair(1).ST) );
            %
            Fc = linmix2( x_dis, obj.lath_solution_pair(1).ST, obj.lath_solution_pair(2).ST );
            disg_delta_hps = delta_hps( Fc );
            disg_rot = axis_angle_rotvec_inclusion( Fc );
            dgo = [x_dis', d_dis, disg_delta_hps, disg_rot];
        end
        
        function glo = get.gl_opt( obj )
            [x_gl,  d_gl] = frob_min_green_lagrange_composite_block( ...
                obj.lath_solution_pair(1).ST, obj.lath_solution_pair(1).ST, obj.shape_vec_opt(1) ); % x_eps(1) );
            %
            Fc = linmix2( x_gl, obj.lath_solution_pair(1).ST, obj.lath_solution_pair(2).ST );
            gl_delta_hps =  delta_hps( Fc );
            rot_svo = axis_angle_rotvec_inclusion( Fc );
            glo = [x_gl', d_gl, gl_delta_hps, rot_svo];
        end 
        
    end % end methods
end % class 

    %% other functions as constructors, get or set functions        

        % deviation of average block habit plane from given crystallographic set/family
        function dhp = delta_hps( Fc )
            [~, ~, ~, ~, h1, h2] = rank_one(Fc, eye(3), 1.e-3, false); % last 'false' is that no lambda_2_warning occurs
            cpps_gamma = all_from_family_perms( [1 1 1] );
            dhp(1) = min_misorientation( cpps_gamma, h1); 
            dhp(2) = min_misorientation( cpps_gamma, h2); 
        end

        % rotation of average (homogenized) deformation of block
        function vec4 = axis_angle_rotvec_inclusion( Fc )
            [~,R] = polardecomposition( Fc );
            vec4 = vrrotmat2vec( R );
            % convert ange to degree
            vec4(4) = rad2deg( vec4(4) );
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








