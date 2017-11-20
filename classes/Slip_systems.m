classdef Slip_systems %< IPS_solution % PET 14.11.17 actually it should be derived but the whole class should be a property of a higher level object!!!
    % slip-solution:
    % class for the used slip-system(s) to modify a stretch tensor to yield an IPS
    % variable amounts of slip systems possible
    
    properties %(Access = public)
        eps_s; % shear magnitudes of normed shear dyads S
        shear_direction; % slip directions of shears (miller indizes)
        slip_normal_plane_vec; % slip plane normals of glide system (miller indizes)
        mirror_plane; % mirror plane of block solution
    end % end of properties
    properties (Dependent)
        stepwidth; % ( planes_between_steps = inverse slip_density (1/g or 1/m) = average nr of slip planes between burgers vector steps 
        max_eps_s;
        min_stepwidth;
    end
    
    methods
        % constructor
        function obj = Slip_systems( varargin ) % F, G, id, eps_0, a, h, Q, LT, g, d1, n1, d2, n2 )
            % All matlab classes have a default constructor with no arguments!
            %
%             if nargin == 0 % =  if isempty(varargin), return;end
%                 % no argument constructor
%                 super_args = {};            
%             end
%             if nargin > 8
%                 super_args = { varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8} };
%             end
%            %
%            obj = obj@IPS_solution( super_args{:} ); % actually only needs: F, G, id, eps_0, d, h, Q, LT
            %
            if nargin > 0 %8
                obj.eps_s = varargin{1}; %{9}; %{1,9};  % initially here was only one 'stepwidth' and it was assummed equal for both slips (see Qi,Khachaturyan 2014 Acta)
                obj.shear_direction = varargin{2}; %{10}; %{1,10};
                obj.slip_normal_plane_vec = varargin{3}; %{11}; %1,11};
            end
            %
            if nargin > 3 %11 % this is only used if a direct averaged block habit plane condition is used as proposed by Qi and Khachaturyan 2014 Acta
                obj.mirror_plane = varargin{4}; %{12}; %{1,11};
            end
        end % end constructors (mostly used to reduce solutions)
        %%        
        function gg = get.stepwidth(obj)
            gg = slip_planes_between_burgerssteps( obj.shear_direction(:,1:3), obj.eps_s, obj.slip_normal_plane_vec(:,1:3), 'cubic'); %TODO generalize to %obj.Bravais_type );
        end
        %
        function max_eps_s = get.max_eps_s(obj)
            max_eps_s = max(obj.eps_s);
        end
        % 
        function min_stepwidth = get.min_stepwidth(obj)
            min_stepwidth = min(obj.stepwidth);
        end
        
    end % end methods
       
end







