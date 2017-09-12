classdef Slip_solution < IPS_solution
    % slip-solution:
    % class for the used slip-system(s) to modify a stretch tensor to yield an IPS
    % variable amounts of slip systems possible
    
    properties (Access = public)
        eps_s = [0.]; % shear magnitudes of normed shear dyads S
        s = [0. 0. 0.]; % slip directions of shears (miller indizes)
        m = [0. 0. 0.]; % slip plane normals of glide system (miller indizes)
        mirror_plane  = [0. 0. 0.]; % mirror plane of block solution
    end % end of properties
    properties (Dependent)
        g; % average nr of slip planes between burgers vector steps 
    end
    
    methods
        % constructor
        function obj = Slip_solution( varargin ) % F, G, id, eps_0, a, h, Q, LT, g, d1, n1, d2, n2 )
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
                obj.eps_s = varargin{1,9};  % initially here was only one g - equal for both slips
                obj.s = varargin{1,10};
                obj.m = varargin{1,11};
            end
            %
            if nargin > 11 % this is only used if a direct averaged block habit plane condition is used as proposed by Qi and Khachaturyan 2014 Acta
                obj.mirror_plane = varargin{1,11};
            end
        end
        
        function gg = get.g(obj)
            gg = slip_planes_between_burgerssteps( obj.s, obj.eps_s, obj.m, 'cubic'); %TODO generalize to %obj.Bravais_type );
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







