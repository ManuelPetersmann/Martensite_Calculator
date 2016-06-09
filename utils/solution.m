classdef solution
  %SOLUTION results of calculation for a slip-system that provides an invariant plane (after kachaturyan)
  %   Detailed explanation goes here
  
  properties (Access = public)
    eps_0 = 0.0 % magnitude of shear
    n1 = [0.0 0.0 0.0] % normal to habit plane (unit vector) - 1st solution
    n2 = [0.0 0.0 0.0] % normal to habit plane (unit vector) - 2nd solution
    a1 = [0.0 0.0 0.0] % shear direction of transformation (unit vector) - 1st solution
    a2 = [0.0 0.0 0.0] % shear direction of transformation (unit vector) - 2nd solution
    Q1 = eye(3) % rotation matrix (for invariant planar match with parent) - 1st solution
    Q2 = eye(3) % rotation matrix (for invariant planar match with parent) - 2nd solution
    g = 0.0 % slip plane spacing (m in paper)
    theta_p % misorientation angle between {111}_gamma and {011}_alpha
  end % end of properties
  
  methods
    % constructor
    function obj = solution(eps_0, a1, a2, n1, n2, Q1, Q2, g, theta_p)
      if nargin > 0
        if nargin >= 8
          obj.eps_0 = eps_0;
          obj.n1 = n1;
          obj.n2 = n2;
          obj.a1 = a1;
          obj.a2 = a2;
          obj.Q1 = Q1;
          obj.Q2 = Q2;
          obj.g = g;
          obj.theta_p = theta_p;
        end
        
        % Ehl: add used slip-systems in solution?
        % composite Bain-block, single-slip system
%         if nargin = 9
%           
%         end
        
        % composite Bain-block, two-slip systems
%         if nargin = 11
%           
%         end        
      end
    end % end of constructor
    
    % function for identifying the solution with lowest misorientation angle 
    % in an array with solutions
    % Ehl: could perhaps be modified to a routine which identifies 
    % the "optimal solution" based on several indicators?!
    function isol_lma = find_lowest_misorientation_angle(solutions)
      
      nsol = size(solutions,2); % get number of solutions in array
      theta_p_min = 999; % init. variable for temp. storage of min. angle
      islo_lma = 0; % init. for index of solution with min. angle
      
      for i = 1 : nsol
        if(solutions(i).theta_p < theta_p_min)
          theta_p_min = solutions(i).theta_p; % update of min. angle
          isol_lma = i; % update of index for solution with min. angle
        end
      end % end of for loop over solutions-array
      
      if(isol_lma == 0) % check, if any reasonable min. angle is found
        error('find_lowest_misorientation_angle: no realistic angle found!')
      end
        
    end % end of function find_lowest_misorientation_angle
    
    % function for formatted output of a solution
    function sol_output(solution)
%       fprintf('==================================================\n')
%       fprintf('= results of calculation for a slip-system that  =\n')
%       fprintf('= provides an invariant plane (after kachaturyan)=\n')
      fprintf('==================================================\n')
      fprintf('magnitude of shear: \n epsilon_0 = %1.4f \n', solution.eps_0)
      fprintf('normal to habit plane (unit vector) - 1st solution: \n n_1 = \n')
      disp (solution.n1)
      fprintf('normal to habit plane (unit vector) - 2nd solution: \n n_2 = \n')
      disp (solution.n2)
      fprintf('shear direction of transformation (unit vector) - 1st solution: \n l_1 = \n')
      disp (solution.a1)
      fprintf('shear direction of transformation (unit vector) - 2nd solution: \n l_2 = \n')
      disp (solution.a2)
      fprintf('rotation matrix (for invariant planar match with parent) - 1st solution: \n R_I1 = \n')
      disp(solution.Q1)
      fprintf('rotation matrix (for invariant planar match with parent) - 2nd solution: \n R_I2 = \n')
      disp(solution.Q2)
      fprintf('slip plane spacing: \n m = %1.4f \n', solution.g)
      fprintf('==================================================\n')
    end % end of function sol_output
    
  end % end of methods
  
end

