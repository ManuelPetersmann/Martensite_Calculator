classdef Solution_array
    % This class provides functions to reduce to and sort for solutions obeying specific
    % criteria e.g. Orientation relations, or shear (slip_density) -, shape strain (lambda_1 - lambda_3) magnitudes
    %  & Martensite - need not to be derived from Martensite class, but is a property of it!
    
    properties
        array;   % entries of array can be objects of type "IPS_solution", "ILS_solution", etc.
        %solutions_available = false; use  isempty(obj.array) instead
        %
        calculation_method; % a string specifying the method to trim the middle eigenvalue to one
        slip_combinations; % nr of possible slip combinations nchoosek, n... total nr of slip systems,
        % k...nr of simultaneously active ones (depending on calculation method)
        selection_criteria; % = containers.Map(); must be initialized int the constructor!!! see https://de.mathworks.com/matlabcentral/answers/331779-using-containers-map-as-a-class-s-property
        % empty containers.map object of "string-criterion" - value pairs c.f. Hashtable, python dict
        cryst_fams; % = containers.Map(); % KS-dirs, NW-dirs, cpps_fam, cp-dir_fam etc...
        % rewrote this class to value class - new variable to store
        % crystallographic families that are added for specific selection criteria
        sorted_after = 'unsorted'; % string specifying criterion array is sorted for
    end
    
    methods
        % constructor for reduction of solutions
        function obj = Solution_array( varargin )
            % All matlab classes have a default constructor with no arguments!
            % This constructor supplied by MATLAB also calls all superclass constructors with no arguments.
            %
            foundnr = 0; % counter for how many matches are found for the construction of constrained solutions

            if nargin > 1 % varargin = { 1- Solution_array object 
                obj.array = varargin{1}; % this first property could even be determined with the code below
%                 ClassName = class( varargin{1}.array(1) );
%                 switch ClassName
%                     case 'IPS_solution'
%                         varargin{1}.array = IPS_solution();
%                     case 'ILS_solution'
%                         varargin{1}.array = ILS_solution();
%                 end
                %% COPY CONSTRUCTOR for solution_array - copies everything execpt - array property
                if isprop(varargin{2},'array')
                    props_to_copy = properties( varargin{2} );
                    for i=1:length( props_to_copy )
                        if ~strcmp(props_to_copy{i},'array')
                            obj.(props_to_copy{i}) = varargin{2}.(props_to_copy{i});
                        end
                    end
                end
            end
            
            %%  ILS specific Stuff
            if nargin == 4 % ILS: 
                for i = 1:size( varargin{2}.array, 2)
                    % abs(lambda2_IPS - 1) > delta_lambda2
                    % rotangle_inclusion   > max_rotation_angle_inclusion
                    if varargin{2}.array(i).(varargin{3}) < varargin{4} % lambda2 or rot_angle_block
                        foundnr = foundnr + 1;
                        obj.array(foundnr) = varargin{2}.array(i);
                    end
                end
            end
            
            %% IPS and partially ILS criteria
            if nargin > 4 % PET. 19.10.17
                %
                if isempty(varargin{2}.selection_criteria)
                    obj.selection_criteria = containers.Map();
                else
                    obj.selection_criteria = varargin{2}.selection_criteria;
                end
                
                % the following creates a new key-value pair if the key is not yet contained,
                % otherwise it overwrites the existing value for the key
                if isempty(varargin{2}.cryst_fams)
                    obj.cryst_fams = containers.Map();
                else
                    obj.cryst_fams = varargin{2}.cryst_fams;
                end
                
                if isempty(varargin{2}.array) % ~varargin{2}.solutions_available
                    error('Empty solutions array given as input for reduction - think about that...');
                end
            end
            
            
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
                for i = 1:size( varargin{2}.array, 2)
                    if strcmp( varargin{5}, 'max')
                        if varargin{2}.array(i).(varargin{3}) < varargin{4} % = maximal tolerated 'eps_ips' value
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    %% the next three are for the slip shear
                    if strcmp( varargin{5}, 'max_shear') 
                        if varargin{2}.array(i).slip.(varargin{3}) < varargin{4} % = maximal tolerated plastic shear value (of normed shear dyads)
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    %
                    if strcmp( varargin{5}, 'shear_sum')
                        if varargin{2}.array(i).slip.(varargin{3}) < varargin{4} % = maximal tolerated plastic shear value (of normed shear dyads)
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    %
                    if strcmp( varargin{5}, 'min')
                        % varargin{2}.array(i).(varargin{3})
                        % varargin{4}
                        % here the condition is only fullfilled if all entries are bigger, i.e. the minimum
                        if varargin{2}.array(i).slip.(varargin{3}) > varargin{4} % = minimal tolerated 'stepwidth' or 'm'/'g' value (inverse '1/m' is called slip density )
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    %%
                    if strcmp( varargin{3}, 'delta_determinant_max')
                        % here varargin{5} is not a string but det(martenstite.U)
                        % determiannt should not change more than some value in varargin{4}
                        if ( abs( det(varargin{2}.array(i).F1) - varargin{5} ) <  varargin{4} )
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                end
            end
            
            
            % to generate minimum misorientation angles with corresponding
            % vector from orientation relations (close packed planes and directions) given as families
            if nargin >= 7 % varargin = { 1 - Type of array entry object, 2 - Solution_array, 3 - characteristic plane or direction-family (mostly CP families),
                % 4 - max tolerated deviation from OR,  5 - dynamic property name for least_misorientation-angle,
                % 6 - dynamic property name for least_misorientation-vector, 7 - dynamic property name for OR family / or 'h' or 'd' for property,
                % 8 - bool for planes = true (otherwise directions), }
                obj.selection_criteria( varargin{5} ) = varargin{4};
                
                for i = 1:size( varargin{2}.array, 2)
                    
                    if isempty(varargin{2}.array(i).added_props)
                        varargin{2}.array(i).added_props = containers.Map();
                    end
                    
                    if nargin == 8
                        obj.cryst_fams( varargin{7} ) = varargin{3};
                        % calculate minimum angle between plane normal
                        % vector (from the given set of vectors) and its
                        % transformed form. (Bain correspondence determines
                        % unmodified correspondence!)
                        % min_misorientation( varargin{3}, varargin{2}.array(i).LT, varargin{8} )
                        [ varargin{2}.array(i).added_props(varargin{5}), varargin{2}.array(i).added_props( varargin{6}) ] = ...
                            min_misorientation( varargin{3}, varargin{2}.array(i).LT, varargin{8} );
                    end
                    
                    
                    if nargin == 7
                        % both if and else do the same just with another syntax!
                        if isprop(varargin{2}.array(i),varargin{7})
                            % misorientation between a crystallographic family
                            % (cp dir or cpp to habit plane or deformation
                            % vector ('h' or 'd' in (7))
                            % add an angle (5) and a vector property -(6) to the
                            % object property 'added_props' if they have not yet been added...
                            %
                            [ varargin{2}.array(i).added_props(varargin{5}), varargin{2}.array(i).added_props(varargin{6}) ] = ...
                                min_misorientation( varargin{3}, varargin{2}.array(i).(varargin{7}) ); % 7 'h', 'd' (IPS_solution property string)
                        else
                            % deviation of prefered invariant line e.g. close packed direction from found habit plane (0 if vector is in the plane)
                            [ varargin{2}.array(i).added_props(varargin{5}), varargin{2}.array(i).added_props(varargin{6}) ] = ...
                                misorientation_vector_and_plane( varargin{3}, varargin{2}.array(i).h );
                        end
                    end
                    %
                    % SORT OUT
                    if varargin{2}.array(i).added_props(varargin{5}) < varargin{4} % varargin{4} = maximal tolerated value
                        foundnr = foundnr + 1;
                        obj.array(foundnr) = varargin{2}.array(i);
                    end
                end
            end
            %
            if nargin > 1
                % after reduction of solutions check if there is at least one non-empty entry in object
                if isempty(obj.array) %(size( obj.array, 2)==1) && isempty(obj.array(1).F1)
                    disp('No Solution fullfilling specified criteria');
                    %    obj.solutions_available = false;
                else
                    disp(['Solutions reduced to : ' , num2str(length(obj.array))] );
                    %    obj.solutions_available = true;
                end
            end
        end
        
        
        %% Sort
        function sorted = sort(obj, prop_name)
            %
            if isempty(obj.array) % ~obj.solutions_available
                error('Empty solutions array cannot be sorted');
            end
            %
            prop_array = zeros(1,size( obj.array ,2) );
            %% TRY WITHOUT FOR LOOP!!! -see Solution_array_composite class
            for i = 1: size( obj.array ,2)
                % value to sort for must be extracted via a loop first
                if isprop(obj.array(1), prop_name) || isprop(obj.array(1).slip, prop_name)
                    if strcmp(prop_name,'stepwidth') % two values - sort for the smaller one 1/stepwidth \propto eps_s (plastic shear magnitude)
                        prop_array(i) = min( obj.array(i).slip.(prop_name) );
                    else
                        prop_array(i) = obj.array(i).(prop_name);
                    end
                else % variables that are not direct properties of solution_array class
                    if isKey(obj.array(i).added_props , prop_name)
                        prop_array(i) = obj.array(i).added_props(prop_name);
                    else
                        error('Solutions cannot be sorted for selection since it is not specified in the selection criteria!');
                    end
                end
                %                [~,idx] = sort( prop_array, 1 );
            end
            %max(prop_array)
            %min(prop_array)
            [~,idx] = sort( prop_array );
            sorted = obj;
            sorted.array = obj.array(idx);
            display(['Solutions sorted in ascending order for property: ' , prop_name ]);
            sorted.sorted_after = prop_name;
        end
        
        
    end % methods
end
     