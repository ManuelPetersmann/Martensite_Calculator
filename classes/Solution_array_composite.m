classdef Solution_array_composite 
    
    properties
        array;   % entries of type "Composite_solution"
        %solutions_available = false;
        %
        calculation_method; % a string specifying the method to trim the middle eigenvalue to one
        %
        selection_criteria; %  = containers.Map(); must be initialized int the constructor!!! see https://de.mathworks.com/matlabcentral/answers/331779-using-containers-map-as-a-class-s-property
        %
        sorted_after = 'unsorted'; % string specifying criterion array is sorted for
    end
    
    methods
        % constructor for reduction of solutions 
        function obj = Solution_array_composite( ) %varargin )
            % All matlab classes have a default constructor with no arguments!
            % This constructor supplied by MATLAB also calls all superclass constructors with no arguments.
            %
            %foundnr = 0; % counter for how many matches are found for the construction of constrained solutions
            obj.array = Composite_solution();
            obj.selection_criteria = containers.Map();
        end
        
        
        %% Sort - copied from Solution_array (laths) maybe adapt sometime
        function [obj,idx]= sort(obj, prop_name)
            %
            if isempty( obj.array )  % ~obj.solutions_available
                error('Empty solutions array cannot be sorted');
            end
            if isprop(obj.array(1),prop_name)
            [~,idx] = sort( obj.array.(prop_name)(3), 1 );
            obj.array = obj.array(idx);
            disp(['Solutions sorted in ascending order for property: ' , prop_name ]);
            obj.sorted_after = prop_name;
            else
                disp('sort string not a property! possibilities are: shape_vec_opt, disg_opt, gl_opt'); 
            end
        end
        
       
    end % methods  
end
     