%% update selection criteria for laths

%handles.martensite.IPS_solutions.array(1).added_props.keys
%handles.martensite.IPS_solutions.cryst_fams
%handles

if handles.martensite.IPS_solutions.solutions_available
    %    
    if handles.asc_number > 0
        updateLog_MartCalc(hObject, handles,'Start reducing solutions after specified criteria. Please wait...');
        
        % criterion 1: Minimum slip plane density
        if(handles.asc_status(1) > 0)
            min_stepwidth = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(1)).Children(2).String);
        end
        
        % Criterion 2: Maximum shape strain
        if(handles.asc_status(2) > 0)
            eps_ips_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(2)).Children(2).String);
        end
        
        % Criterion 3: Maximum misorientation of CPPs {110}_alpha and {111}_gamma
        if(handles.asc_status(3) > 0)
            theta_CPPs_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(3)).Children(2).String);
        end
        
        % Criterion 4: Maximum misorientation of block HP to {111}_gamma
        % note 557 is 9.4° from 111 ! therefore this high tolerance!
        % Angle between 111 and 557 habit plane
        % acos( dot([1. 1. 1.], [5. 5. 7.])/(sqrt(3)*sqrt(99) ) ) = 9.4 degree
        if(handles.asc_status(4) > 0)
            theta_h_to_cpp = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(4)).Children(2).String);
        end
        
        % Criterion 5: Maximum deviation from KS OR
        % The peak of OR distribution is normally between KS and NW and
        % these two are 5.25 apart - hence these tolerances
        % 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
        if(handles.asc_status(5) > 0)
            theta_KS_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(5)).Children(2).String);
        end
        
        % Criterion 6 has been chosen: Maximum deviation from NW OR
        %'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
        if(handles.asc_status(6) > 0)
            theta_NW_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(6)).Children(2).String);
        end
        
        % PET: 12.10.17
        % Criterion 8 - deviation of preferred ILS direction from invariant plane (0 if vector is in plane)
        if(handles.asc_status(7) > 0)
            theta_max_ILSdir_to_h = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(7)).Children(2).String);
        end
        
        %    % PET: Moved this criterion from 5 to 8 and commented it! use tolerance
        %    % of 0.001 per default! everything other is not reasonable...
        %    % Criterion 8: Maximum deviation of determinant det(F) of transformation
        %     if(handles.asc_status(8) > 0)
        %         delta_determinant_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(5)).Children(2).String);
        %     end
        % Also removed in GUIDE - formerly - Maximum deviation of theoretical volume change of Bain strain
        
        %% reduce solutions
        %
        red_sols = handles.martensite.IPS_solutions;
        for criterion = 1:handles.asc_number
            if ~red_sols.solutions_available
                break
            end
            %
            % here conditional assignments are used c.f. c++ - (exp1) ? exp2 : exp3  - If exp1 == true, then exp2 else exp3
            % to check wheter the no solutions are found (selection criteria too strict)
            switch handles.asc_list( criterion )
                case 1 % stepwidth
                    red_sols = Solution_array( Slip_solution, red_sols, 'stepwidth', min_stepwidth, 'min');
                    crit = [' for a stepwidth > ',num2str(min_stepwidth)];
                case 2 % eps_ips
                    red_sols =   Solution_array( Slip_solution, red_sols, 'eps_ips', eps_ips_max, 'max' );
                    crit = [' for a shape strain  < ',num2str(eps_ips_max)];
                case 3 % theta_CPPs
                    red_sols =    Solution_array( Slip_solution, red_sols, handles.austenite.CPPs, theta_CPPs_max, 'theta_CPPs', 'closest_CPPs', 'cpps_gamma', true);
                    crit = [' for deviation from ideal CP relation  < ',num2str(theta_CPPs_max),'°'];
                case 4 % theta_h_to_cpp
                    red_sols =    Solution_array( Slip_solution, red_sols, handles.austenite.CPPs, theta_h_to_cpp, 'theta_h_to_CPP', 'closest_h_to_CPP', 'h');
                    crit = [' for a habit plane misorientation to {111}_aust  < ',num2str(theta_h_to_cpp),'°'];
                case 5
                    red_sols =    Solution_array( Slip_solution, red_sols, handles.austenite.CP_dirs, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );
                    crit = [' for a maximum deviation angle between ideal KS-direction-parallelism  < ',num2str(theta_KS_max),'°'];
                case 6
                    red_sols =    Solution_array( Slip_solution, red_sols, handles.NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);
                    crit = [' for a maximum deviation angle between ideal NW-direction-parallelism  < ',num2str(theta_NW_max),'°'];
                case 7 % theta_max_ILSdir_to_h
                    red_sols =   Solution_array( Slip_solution, red_sols, handles.austenite.CP_dirs, theta_max_ILSdir_to_h, 'theta_preferred_ILSdir_to_h', 'closest_ILSdir_to_h','KS');
                    crit = [' for a maximum deviation angle of preferred invariant line from invariant habit plane < ',num2str(theta_max_ILSdir_to_h)];
               % case 8
                    % red_sols = Solution_array( Slip_solution, red_sols, 'delta_determinant_max', delta_determinant_max,  det(handles.martensite.U));
                    %  crit = [' for (non-physical) volume change  > ',num2str(delta_determinant_max)];
            end
            %
            if red_sols.solutions_available
                updateLog_MartCalc(hObject, handles,['Solutions reduced to : ', num2str(length(red_sols.array)), crit ] );
                handles.reduced_solutions = red_sols;
            else
                updateLog_MartCalc(hObject, handles,'No Solution fullfilling all specified criteria. Solution reduction stopped before next active selection criterion (asc). See asc list.');
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
        updateLog_MartCalc(hObject, handles, 'Note: No filters for solutions specified.');
    end
else
    updateLog_MartCalc(hObject, handles, 'No lath solutions available.');
end % end if .solutions_available = true

