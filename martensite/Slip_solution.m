classdef Slip_solution < IPS_solution
    % slip-solution:
    % class for the found slip-system(s), slip amount - g solutions, for the
    % modification of the stretch tensor to yield an IPS
    % maximal two slip systems possible
    
    properties (Access = public)
        g = 0.0 % slip plane spacing (m in paper)
        d1 = [0. 0. 0.];
        d2 = [0. 0. 0.];
        n1 = [0. 0. 0.];
        n2 = [0. 0. 0.];        
    end % end of properties
    
    methods
        % constructor
        function obj = Slip_solution( varargin ) % F, G, id, eps_0, a, h, Q, LT, g, d1, n1, d2, n2 )
            %
            if nargin == 0 % =  if isempty(varargin), return;end
                % no argument constructor
                super_args = {};            
            end
            if nargin > 9
                super_args = { varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8} };
            end
            %
            obj = obj@IPS_solution( super_args{:} ); % actually only needs: F, G, id, eps_0, a, n, Q, LT
            %
            if nargin > 9
                obj.g = varargin{1,9};
                obj.d1 = varargin{1,10};
                obj.n1 = varargin{1,11};
            end
            %
            if nargin > 11
                obj.d2 = varargin{1,12};
                obj.n2 = varargin{1,13};
            end
        end
        
        % function for formatted output of a solution
%         function sol_output(solution)
%             %       fprintf('==================================================\n')
%             %       fprintf('= results of calculation for a slip-system that  =\n')
%             %       fprintf('= provides an invariant plane (after kachaturyan)=\n')
%             
%             format long
%             
%             fileID = fopen('Habitplane_evaluation.txt', 'a');
%             
%             format_m = 'm   = (%d, %d, %d) \t';
%             format_l = 'l   = [%d, %d, %d] \n';
%             format_g = 'g = %7.4f \t';
%             
%             fprintf(fileID,'shear system nr: %d \n', idx(i,1) );
%             fprintf(fileID,format_m, mm);
%             fprintf(fileID,format_l, ll);
%             fprintf(fileID,format_g, g);
%             
%             fprintf('slip plane spacing: \n m = %1.4f \n \n', solution.g)
%             
%         end
        
    end % end methods
       
end







