classdef Solution_array < dynamicprops % subclass of handle class
    % This class provides functions to sort for solutions obeying specific
    % criteria e.g. Orientation relations, or shear (slip_density) -, shape strain (lambda_1 - lambda_3) magnitudes
    %  & Martensite - need not to be derived from Martensite class, but is
    %  a property of it!
    % must be derived from dynamicprops because subclass is derived from it
    properties
        array;   % entries of array can be objects of type "IPS_solution", "Slip_solution", etc.
        no_solutions_available = false;
        %
        calculation_method; % a string specifying the method to trim the middle eigenvalue to one
        slip_combinations; % nr of possible slip combinations nchoosek, n... total nr of slip systems, 
        % k...nr of simultaneously active ones (depending on calculation method)
        selection_criteria = containers.Map; % containers.map object of "string-criterion" - value pairs c.f. Hashtable, python dict
        sorted_after; % string specifying criterion array is sorted for
    end
    
    methods
        % constructor
        function obj = Solution_array( varargin )
            % All matlab classes have a default constructor with no arguments!
            foundnr = 0; % counter for how many matches are found for the construction of constrained solutions
            if nargin == 1 % varargin = { 1-Type of array entry object }
                % initalize array type e.g. IPS_solutions -or-> Slip_solution
                obj.array = varargin{1}; % here varargin{1} should be the object_type of class property array
            end
            %
            if (nargin > 1) && varargin{2}.no_solutions_available
                error('No solutions fullfilling specified selection criteria.')
            end
            %
            if nargin == 5 % to construct subarrays with minimal/maximal 'slip_density' and 'eps_ips' values
                % as well as ones with specified maximum change of determinant
                % varargin = {1-Type of array entry object, 
                % 2 - Solution_array object, 3 - object property, 4 -extremum
                % value for reduction of solutions, 5 - string specifying
                % if value from 4 is the allowed maximum 'max or minimum 'min'
                %
                %
                % the following creates a new key-value pair if the key is
                % not yet contained, otherwise it overwrites the existing
                % value for the key
                obj.selection_criteria( varargin{3} ) = varargin{4};
                %
                obj.array = varargin{1};
                for i = 1:size( varargin{2}.array, 2)
                    if strcmp( varargin{5}, 'max')
                        if varargin{2}.array(i).(varargin{3}) < varargin{4} % = maximal tolerated 'eps_ips' value
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    if strcmp( varargin{5}, 'min')
                        if varargin{2}.array(i).(varargin{3}) > varargin{4} % = minimal toleratred 'slip_density' value
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    if strcmp( varargin{3}, 'delta_determinant_max')
                        % determiannt should not change more than some value in varargin{4}
                        if ( abs( det(varargin{2}.array(i).F1) - varargin{5} ) <  varargin{4} )
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                end
            end
            %
            % to generate minimum misorientation angles with corresponding
            % vector from orientation relations (close packed planes and directions) given as families
            if nargin >= 6 % varargin = { 1 - Type of array entry object, 2 - Solution_array, 3 - characteristic plane or direction-family (mostly CP families), 
                % 4 - max tolerated deviation from OR,  5 - dynamic property name for least_misorientation-angle, 
                % 6 - dynamic property name for least_misorientation-vector, 7 - dynamic property name for OR family / or 'h' or 'd' for property,  
                %8- bool for planes = true (otherwise directions), }
                obj.array   = varargin{1};  
                %
                % the following creates a new key-value pair if the key is
                % not yet contained, otherwise it overwrites the existing
                % value for the key
                obj.selection_criteria( varargin{5} ) = varargin{4};
                %
                if nargin == 8
                    if ~isprop(obj,varargin{7}) % if property not already added add it dynamically
                        obj.addprop( varargin{7} );
                    end
                    obj.( varargin{7} ) = varargin{3};
                end
                for i = 1:size( varargin{2}.array, 2)
                    % add an angle (5) and a vector property -(6) to the object
                    % dynamically if they have not yet been added...
                    if ~isprop(varargin{2}.array(i),varargin{5}) 
                        varargin{2}.array(i).addprop( varargin{5} );
                    end
                    if ~isprop(varargin{2}.array(i),varargin{6}) 
                        varargin{2}.array(i).addprop( varargin{6} );
                    end
                    %
                    if nargin == 8
                        % calculate minimum angle between plane normal
                        % vector (from the given set of vectors) and its
                        % transformed form. (Bain correspondence determines
                        % unmodified correspondence!)
                        [ varargin{2}.array(i).(varargin{5}), varargin{2}.array(i).(varargin{6}) ] = ...
                            min_misorientation( varargin{3}, varargin{2}.array(i).LT, varargin{8} );
                    end
                    %
                    if nargin == 7
                        % misorientation between a crystallographic family
                        % (cp dir or cpp to habit plane or deformation
                        % vector ('h' or 'd' in (7))
                        [ varargin{2}.array(i).(varargin{5}), varargin{2}.array(i).(varargin{6}) ] = ...
                            min_misorientation( varargin{3}, varargin{2}.array(i).(varargin{7}) ); % 7 'h', 'd' or 'e1'
                    end
                    %
                    if nargin == 6
                        % deviation of prefered invariant line e.g. close packed direction from found habit plane (0 if vector is in the plane) 
                        [ varargin{2}.array(i).(varargin{5}), varargin{2}.array(i).(varargin{6}) ] = ...
                            misorientation_vector_and_plane( varargin{3}, varargin{2}.array(i).h ); 
                    end
                    % SORT OUT 
                    if varargin{2}.array(i).(varargin{5}) < varargin{4} % varargin{4} = maximal tolerated value
                        foundnr = foundnr + 1;
                        obj.array(foundnr) = varargin{2}.array(i);
                    end                   
                end
            end
            %
            if nargin > 1
                % after reduction of solutions check if there is at least one non-empty entry in object
                if (size( obj.array, 2)==1) && isempty(obj.array(1).F1)
                    disp('No Solution fullfilling specified criteria');
                    obj.no_solutions_available = true;
                else
                    disp(['Solutions reduced to : ' , num2str(length(obj.array)) ]);
                end
            end
        end

            
        %% Sort 
        function [obj,idx]= sort(obj, prop_name)
            if isprop(obj, prop_name)
                [~,idx] = sort( obj.array.(prop_name)  ); % key here is the obj.('prop_name') - 
                % see http://de.mathworks.com/help/matlab/matlab_prog/generate-field-names-from-variables.html
            else % dynamic properties not directly accesible like above therefore the array that is sorted 
                 % must be extracted via a loop first
                prop_array = zeros(1,size( obj.array ,2) );
                for i = 1: size( obj.array ,2)
                    prop_array(i) = obj.array(i).(prop_name);
                end
                [~,idx] = sort( prop_array );
            end
            obj.array = obj.array(idx);
            display(['Solutions sorted in ascending order for property: ' , prop_name ]);
            obj.sorted_after = prop_name;
        end
        
       
    end % methods  
end
     