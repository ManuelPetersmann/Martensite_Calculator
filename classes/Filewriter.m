        function write_input_parameters(obj, martobj, austobj)
                        %format long
                        FileID = fopen(filename,'w');
            fprintf(obj.FileID, '%s\n', text_str);
            % write date
            datestr(now,'HH:MM:SS.FFF')
            % write lattice parameters
            austobj.Lp(1)
            martobj.Lp(1)
            % write slip plane/direction families of phases
            martobj.considered_plasticity
            martobj.slip_planes;
            martobj.slip_directions;
            austobj.slip_planes;
            austobj.slip_directions;
            % write correspondence matrix (coordinate trans)
            
            % write Bases
            handles.martensite.my_base
            handles.austenite.my_base
            
            fclose(FileID);
        end
        
        %---------------------------------------------------------------
        function write_calc_specs(obj,)
            % 5
            obj.calculation_method
            

            

            
            fileID = fopen('Habitplane_evaluation.txt', 'a');
            
            format_m = 'm   = (%d, %d, %d) \t';
            format_l = 'l   = [%d, %d, %d] \n';
            format_g = 'g = %7.4f \t';
            
            fprintf(fileID,'shear system nr: %d \n', idx(i,1) );
            fprintf(fileID,format_m, mm);
            fprintf(fileID,format_l, ll);
            fprintf(fileID,format_g, g);
            
            fprintf('slip plane spacing: \n m = %1.4f \n \n', solution.g)
            % format_variant = ' U%d = %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n \n';
            % fprintf(fileID,format_variant,variantnr,B);
            %
            % format_e = 'eigs = %7.4f %7.4f %7.4f \n';
            % format_mod_variant = ' U_mod = %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n     %7.4f %7.4f %7.4f \n \n';
            % fprintf(fileID,format_mod_variant, mod_variant);
        end



       
        
