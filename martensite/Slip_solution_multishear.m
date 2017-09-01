classdef Slip_solution_multishear < IPS_solution
    % slip-solution:
    % class for the found slip-system(s), slip amount - g solutions, for the
    % modification of the stretch tensor to yield an IPS
    % maximal two slip systems possible
    
    properties (Access = public)
        shear_mags = []; % slip plane spacings of shears
        ds = []; % shear directions
        ns = []; % according shear plane normals    
    end % end of properties
    
    methods
        % constructor
        function obj = Slip_solution_multishear( varargin ) % F, G, id, eps_0, a, h, Q, LT,   ---  g, d1, n1, d2, n2 )
            %
            if nargin == 0 % =  if isempty(varargin), return;end
                % no argument constructor
                super_args = {};            
            end
            if nargin > 8
                super_args = { varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8} };
            end
            %
            obj = obj@IPS_solution( super_args{:} ); % actually only needs: F, G, id, eps_0, a, n, Q, LT
            %
            if nargin > 8
                obj.shear_mags = varargin{1,9};  
                obj.ds = varargin{1,10};
                obj.ns = varargin{1,11};
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







