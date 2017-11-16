%% update selection criteria for BLOCKS

if ~isempty(handles.martensite.IPS_solutions.array) || %solutions_available
    %    
    if handles.asc_number_blocks > 0
        
        det_tol = num2str(handles.handles.edit_minors_det.String);
        if( det_tol <= 1.e-3 )
            XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX = str2num(handles.edit_minors_det.String);
        else
            updateLog_MartCalc(hObject, handles,'A higher tolerance than 1.e-3 for minors relations is not allowed.');
            handles.input_status = false;
        end
        
        cof_tol = num2str(handles.handles.edit_minors_cof.String);
        if( cof_tol <= 1.e-3 )
            XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX = str2num(handles.edit_minors_cof.String);
        else
            updateLog_MartCalc(hObject, handles,'A higher tolerance than 1.e-3 for minors relations is not allowed.');
            handles.input_status = false;
        end
        
        
        
        
 %% TODO  CONTINUE
 
 
        % criterion 1: Minimum slip plane density
        if(handles.asc_status_IPS(1) > 0)
            min_stepwidth = str2num(handles.pan_asc_IPS.Children(size(handles.pan_asc_IPS.Children,1)+1-handles.asc_status_IPS(1)).Children(2).String);
        end
        
        % Criterion 2: Maximum shape strain
        if(handles.asc_status_IPS(2) > 0)
            eps_ips_max = str2num( handles.pan_asc_IPS.Children( size(handles.pan_asc_IPS.Children,1)+1-handles.asc_status_IPS(2) ).Children(2).String );
        end
        
        % Criterion 3: Maximum misorientation of CPPs {110}_alpha and {111}_gamma
        if(handles.asc_status_IPS(3) > 0)
            theta_CPPs_max = str2num(handles.pan_asc_IPS.Children(size(handles.pan_asc_IPS.Children,1)+1-handles.asc_status_IPS(3)).Children(2).String);
        end
        
        %% reduce solutions
        updateLog_MartCalc(hObject, handles,'Start reducing solutions after specified criteria. Please wait...');
        %
        red_sols = handles.martensite.IPS_solutions;
        for criterion = 1:handles.asc_number_IPS
            if isempty(red_sols.array) % ~red_sols.solutions_available
                break
            end
            %
            % here conditional assignments are used c.f. c++ - (exp1) ? exp2 : exp3  - If exp1 == true, then exp2 else exp3
            % to check wheter the no solutions are found (selection criteria too strict)
            switch handles.asc_list_IPS( criterion )
                case 1 % stepwidth
                    red_sols = Solution_array( IPS_solution, red_sols, 'stepwidth', min_stepwidth, 'min');
                    crit = [' for a stepwidth > ',num2str(min_stepwidth)];
                case 2 % eps_ips
                    red_sols =   Solution_array( IPS_solution, red_sols, 'eps_ips', eps_ips_max, 'max' );
                    crit = [' for a shape strain  < ',num2str(eps_ips_max)];
                case 3 % theta_CPPs
                    red_sols =    Solution_array( IPS_solution, red_sols, handles.austenite.CPPs, theta_CPPs_max, 'theta_CPPs', 'closest_CPPs', 'cpps_gamma', true);
                    crit = [' for deviation from ideal CP relation  < ',num2str(theta_CPPs_max),'°'];
                case 4 % theta_h_to_cpp
                    red_sols =    Solution_array( IPS_solution, red_sols, handles.austenite.CPPs, theta_h_to_cpp, 'theta_h_to_CPP', 'closest_h_to_CPP', 'h');
                    crit = [' for a habit plane misorientation to {111}_aust  < ',num2str(theta_h_to_cpp),'°'];
                case 5
                    red_sols =    Solution_array( IPS_solution, red_sols, handles.austenite.CP_dirs, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );
                    crit = [' for a maximum deviation angle between ideal KS-direction-parallelism  < ',num2str(theta_KS_max),'°'];
                case 6
                    red_sols =    Solution_array( IPS_solution, red_sols, handles.NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);
                    crit = [' for a maximum deviation angle between ideal NW-direction-parallelism  < ',num2str(theta_NW_max),'°'];
                case 7 % theta_max_ILSdir_to_h
                    red_sols =   Solution_array( IPS_solution, red_sols, handles.austenite.CP_dirs, theta_max_ILSdir_to_h, 'theta_preferred_ILSdir_to_h', 'closest_ILSdir_to_h','KS');
                    crit = [' for a maximum deviation angle of preferred invariant line from invariant habit plane < ',num2str(theta_max_ILSdir_to_h)];
               % case 8
                    % red_sols = Solution_array( IPS_solution, red_sols, 'delta_determinant_max', delta_determinant_max,  det(handles.martensite.U));
                    %  crit = [' for (non-physical) volume change  > ',num2str(delta_determinant_max)];
            end
            %
            if isempty(red_sols.array) %red_sols.solutions_available
                updateLog_MartCalc(hObject, handles,'No Solution fullfilling all specified criteria. Solution reduction stopped before next active selection criterion (asc). See asc list.');
            else
                updateLog_MartCalc(hObject, handles,['Solutions reduced to : ', num2str(length(red_sols.array)), crit ] );
                handles.reduced_solutions = red_sols;
            end
            %
            guidata(hObject, handles);
            %handles.reduced_solutions.cryst_fams.keys
            
        end % end for       
%         red_sols.array(1).F1
%         ~isempty(red_sols.array(1).F1)
%         ~(size( red_sols.array, 2)==1)
%         handles.reduced_solutions       
%         handles.reduced_solutions.selection_criteria.keys
%         handles.reduced_solutions.array(1).added_props.keys      
%         handles.reduced_solutions.array(1).F1
%         handles.reduced_solutions.array(2).F1
        
        updateLog_MartCalc(hObject, handles, 'Filtering of IPS solutions after specified criteria completed.');
    else
        updateLog_MartCalc(hObject, handles, 'Note: No filters for IPS solutions specified.');
    end
else
    updateLog_MartCalc(hObject, handles, 'No IPS lath solutions available.');
end % end if isempty -- .solutions_available = true

