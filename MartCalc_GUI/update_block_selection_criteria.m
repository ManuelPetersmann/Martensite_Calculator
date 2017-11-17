%% update selection criteria for BLOCKS


if ( ( ~isempty(handles.martensite.IPS_solutions.array) ) && (handles.popup_calc_lath_level.Value == 1) )  || ...
   ( ( ~isempty(handles.martensite.ILS_solutions.array) ) && (handles.popup_calc_lath_level.Value == 2) )
    
    % Minors reltation tolerances
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
    
    if handles.asc_number_blocks > 0
        
        % criterion 1: rotation angle of (average) block deformation
        if(handles.asc_status_blocks(1) > 0)
            max_rot_angle_block = str2num(handles.pan_asc_blocks.Children(size(handles.pan_asc_blocks.Children,1)+1-handles.asc_status_blocks(1)).Children(2).String);
        end
        
        % Criterion 2: Deviation of IPS condition (lambda2 from 1)
        if(handles.asc_status_blocks(2) > 0)
            lambda2_tol_block_aust = str2num( handles.pan_asc_blocks.Children( size(handles.pan_asc_blocks.Children,1)+1-handles.asc_status_blocks(2) ).Children(2).String );
        end
        
        % Criterion 3: If crit 2 is valid - defiation of block HP to {111}_aust
        if(handles.asc_status_blocks(3) > 0)
            block_hp_cp_aust_tol = str2num(handles.pan_asc_blocks.Children(size(handles.pan_asc_blocks.Children,1)+1-handles.asc_status_blocks(3)).Children(2).String);
        end
        
        %% reduce solutions
        updateLog_MartCalc(hObject, handles,'Start reducing solutions after specified criteria. Please wait...');
        %
        red_sols = handles.martensite.IPS_solutions;
        for criterion = 1:handles.asc_number_blocks
            if isempty(red_sols.array) % ~red_sols.solutions_available
                break
            end
            %
            % here conditional assignments are used c.f. c++ - (exp1) ? exp2 : exp3  - If exp1 == true, then exp2 else exp3
            % to check wheter the no solutions are found (selection criteria too strict)
            switch handles.asc_list_blocks( criterion )
                case 1 % stepwidth
                    red_sols = Solution_array( some_solution, red_sols, 'angle_rotvec_inclusion', max_rot_angle_block');
                    crit = [' for a stepwidth > ',num2str(min_stepwidth)];
                case 2 % eps_ips
                    red_sols =   Solution_array( some_solution, red_sols, 'eps_ips', lambda2_tol_block_aust);
                    crit = [' for a shape strain  < ',num2str(eps_ips_max)];
                case 3 % theta_h_to_cpp
                    red_sols =    Solution_array( IPS_solution, red_sols, handles.austenite.CPPs, theta_h_to_cpp, 'theta_h_to_CPP', 'closest_h_to_CPP', 'h');
                    crit = [' for a habit plane misorientation to {111}_aust  < ',num2str(theta_h_to_cpp),'°'];
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
            
        end % end for
        
        updateLog_MartCalc(hObject, handles, 'Filtering of Block solution build from IPS lath solutions after specified criteria completed.');
    else
        updateLog_MartCalc(hObject, handles, 'Note: No filters for Block solution build from IPS lath solutions specified.');
    end
else
    updateLog_MartCalc(hObject, handles, 'No IPS lath solutions available.');
end 

