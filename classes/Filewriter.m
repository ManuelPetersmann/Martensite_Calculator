classdef Filewriter < handle
    % Basic class for writing see matlab_oop.pdf - 3-20
    properties (Access = private)
        FileID
    end
    
    methods
        %constructor
        function obj = Filewriter(filename)
            obj.FileID = fopen(filename,'w');
        end
        
        function writeToFile(obj,text_str)
            fprintf(obj.FileID, '%s\n', text_str);
        end
        
        function delete(obj)
            fclose(obj.FileID);
        end
    end
end