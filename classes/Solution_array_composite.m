classdef Solution_array_composite 
    
    properties
        array;   % entries of array can be objects of type "IPS_solution", "Slip_solution", etc.
        %solutions_available = false;
        %
        calculation_method; % a string specifying the method to trim the middle eigenvalue to one
        mixing_tolerances; %  = containers.Map(); must be initialized int the constructor!!! see https://de.mathworks.com/matlabcentral/answers/331779-using-containers-map-as-a-class-s-property
        % TODO sorted_after = 'unsorted'; % string specifying criterion array is sorted for
    end
    
    methods
        % constructor for reduction of solutions 
        function obj = Solution_array_composite( ) %varargin )
            % All matlab classes have a default constructor with no arguments!
            % This constructor supplied by MATLAB also calls all superclass constructors with no arguments.
            %
            %foundnr = 0; % counter for how many matches are found for the construction of constrained solutions
            obj.array = Composite_solution();
            obj.mixing_tolerances = containers.Map();
        end
        
        
            
        %% Sort - copied from Solution_array (laths) maybe adapt sometime
%         function [obj,idx]= sort(obj, prop_name)
%             %
%             if ~obj.solutions_available
%                 error('Empty solutions array cannot be sorted');
%             end
%             %
%             prop_array = zeros(1,size( obj.array ,2) );
%             for i = 1: size( obj.array ,2)
%                 % value to sort for must be extracted via a loop first
%                 if isprop(obj.array, prop_name)
%                     if strcmp(prop_name,'stepwidth') % two values - sort for the smaller one 1/stepwidth \propto eps_s (plastic shear magnitude)
%                         prop_array(i) = min(obj.array(i).(prop_name));
%                     else
%                         prop_array(i) = obj.array(i).(prop_name);
%                     end
%                 else % variables that are not direct properties of solution_array class
%                     if isKey(obj.array(i).added_props , prop_name)
%                         prop_array(i) = obj.array(i).added_props(prop_name);
%                     else
%                         error('Solutions cannot be sorted for selection since it is not specified in the selection criteria!');
%                     end
%                 end
%                 [~,idx] = sort( prop_array, 1 );
%             end
%             obj.array = obj.array(idx);
%             display(['Solutions sorted in ascending order for property: ' , prop_name ]);
%             obj.sorted_after = prop_name;
%         end
        
       
    end % methods  
end
     