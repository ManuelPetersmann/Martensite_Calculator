classdef Solution_array < dynamicprops % subclass of handle class
    % This class provides functions to sort for solutions obeying specific
    % criteria e.g. Orientation relations, or shear (g) -, shape strain (lambda_1 - lambda_3) magnitudes
    % & Martensite - need not to be derived from Martensite class
    % must be derived from dynamicprops because subclass is derived from it
    properties
        array;   % entries of array can be objects of type "IPS_solution", "Slip_solution", etc.
    end
    
    methods
        % constructor
        function obj = Solution_array( varargin )
            foundnr = 0; % counter for how many matches are found for the construction of constrained solutions
            if nargin == 1 % varargin = { 1-Type of array entry object }
                % initalize array type
                obj.array = varargin{1}; % here varargin{1} should be the object_type of array(i)
            end
            %
            if nargin == 5 % to construct subarrays with minimal/maximal g and eps values
                % as well as ones with specified maximum change of determinant
                % varargin = {1-Type of array entry object, 2
                % -Solution_array object, 3- object property, 4-extremum
                % value for reduction of solutions, 5- string specifying
                % if value from 4 is the allowed maximum 'max or minimum 'min'
                obj.array = varargin{1};
                for i = 1:size( varargin{2}.array, 2)
                    if strcmp( varargin{5}, 'max')
                        if varargin{2}.array(i).(varargin{3}) < varargin{4} % = maximal toleratred eps value
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    if strcmp( varargin{5}, 'min')
                        if varargin{2}.array(i).(varargin{3}) > varargin{4} % = minimal toleratred g value
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                    if strcmp( varargin{3}, 'det')
                        % determiannt should not change more than some value in varargin{4}
                        if ( det(varargin{2}.array(i).F_ips) - varargin{5}) <  varargin{4}
                            foundnr = foundnr + 1;
                            obj.array(foundnr) = varargin{2}.array(i);
                        end
                    end
                end
            end
            %
            % to generate minimum misorientation angles with corresponding
            % vector from orientation relations given as families
            if nargin >= 7 % varargin = { 1- Type of array entry object, 2- Solution_array, 3- OR-family, 4-max tolerated deviation from OR, 
                obj.array   = varargin{1}; % 5- dynamic property name for least_misorientation-angle, 6- dynamic property name for least_misorientation-vector, 
                if nargin == 8             % 7- dynamic property name for OR family / or 'h' or 'a' for property, 8- bool for planes = true (otherwise directions), }
                        obj.addprop( varargin{7} );
                        obj.( varargin{7} ) = varargin{3};
                end
                for i = 1:size( varargin{2}.array, 2)
                    varargin{2}.array(i).addprop( varargin{5} );
                    varargin{2}.array(i).addprop( varargin{6} );
                    if nargin == 8
                        [ varargin{2}.array(i).(varargin{5}), varargin{2}.array(i).(varargin{6}) ] = min_misorientation( varargin{3}, varargin{2}.array(i).LT, varargin{8} );
                    end
                    if nargin == 7
                        % closest austenite direction or plane to 'n' or 'a'
                        [ varargin{2}.array(i).(varargin{5}), varargin{2}.array(i).(varargin{6}) ] = min_misorientation( varargin{3}, varargin{2}.array(i).(varargin{7}) ); 
                    end
                    %
                    if varargin{2}.array(i).(varargin{5}) < varargin{4} % varargin{4} = maximal tolerated value
                        foundnr = foundnr + 1;
                        obj.array(foundnr) = varargin{2}.array(i);
                    end                   
                end
            end
            %
            if nargin > 1
             display(['  current solutions: ' , num2str(length(obj.array)) ]);
            end
        end

            
        %% Sort 
        function [obj,idx]= sort(obj, prop_name)
            if isprop(obj, prop_name)
                [~,idx] = sort( obj.array.(prop_name)  ); % key here is the obj.('prop_name') - see http://de.mathworks.com/help/matlab/matlab_prog/generate-field-names-from-variables.html
            else % dynamic properties not directly accesible like above therefore the array that is sorted must be extracted via a loop first
                prop_array = zeros(1,size( obj.array ,2) );
                for i = 1: size( obj.array ,2)
                    prop_array(i) = obj.array(i).(prop_name);
                end
                [~,idx] = sort( prop_array );
            end
            obj.array = obj.array(idx);
        end
        
       
    end % methods  
end
     