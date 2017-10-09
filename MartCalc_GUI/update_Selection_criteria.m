%% update selection criteria for laths

if handles.asc_number > 0
    updateLog_MartCalc(hObject, handles,'Start reducing solutions after specified criteria.');
    
    % criterion 1: Minimum slip plane density
    if(handles.asc_status(1) > 0)
        g_min = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(1)).Children(2).String);
    end
    
    % Criterion 2: Maximum shape strain
    if(handles.asc_status(2) > 0)
        eps_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(2)).Children(2).String);
    end
    
    % Criterion 3: Maximum misorientation of CPPs {110}_alpha and {111}_gamma
    if(handles.asc_status(3) > 0)
        theta_CPP_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(3)).Children(2).String);
    end
    
    % Criterion 4: Maximum misorientation of block HP to {111}_gamma
    % note 557 is 9.4° from 111 ! therefore this high tolerance!
    % Angle between 111 and 557 habit plane
    % acos( dot([1. 1. 1.], [5. 5. 7.])/(sqrt(3)*sqrt(99) ) ) = 9.4 degree
    if(handles.asc_status(4) > 0)
        theta_n_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(4)).Children(2).String);
    end
    
    % Criterion 5: Maximum deviation of determinant det(F) of transformation
    if(handles.asc_status(5) > 0)
        delta_determinant_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(5)).Children(2).String);
    end
    
    % Criterion 6: Maximum deviation from KS OR
    % The peak of OR distribution is normally between KS and NW and
    % these two are 5.25 apart - hence these tolerances
    % 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
    if(handles.asc_status(6) > 0)
        theta_KS_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(6)).Children(2).String);
    end
    
    % Criterion 7 has been chosen: Maximum deviation from NW OR
    %'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
    if(handles.asc_status(7) > 0)
        theta_NW_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(7)).Children(2).String);
    end
    
    
    %% reduce solutions
    reduced_solutions = handles.martensite.IPS_solutions;
    %
    for criterion = 1:handles.asc_number
        if (size( reduced_solutions.array, 2)==1) && isempty(reduced_solutions.array(1).F1)
            updateLog_MartCalc(hObject, handles,'No Solution fullfilling specified criteria');
            reduced_solutions.no_solutions_available = true;
            break
        end
        %
        switch handles.asc_list( criterion )
            case 1
                reduced_solutions = Solution_array( Slip_solution(), reduced_solutions, 'slip_density', g_min, 'min');
                updateLog_MartCalc(hObject, handles, ['Solutions reduced to : ', num2str(length(reduced_solutions.array)) ' for a slip density > ',num2str(g_min)] );
            case 2
                reduced_solutions = Solution_array( Slip_solution(), reduced_solutions, 'eps_ips', eps_max, 'max' );
                updateLog_MartCalc(hObject, handles, ['Solutions reduced to : ', num2str(length(reduced_solutions.array)) ' for a shape strain  < ',num2str(eps_max)] );
            case 3
                reduced_solutions = Solution_array( Slip_solution(), reduced_solutions, handles.cpps_gamma, theta_CPP_max, 'theta_CPP', 'closest_to_cpp', 'cpps_gamma', true);
                updateLog_MartCalc(hObject, handles, ['Solutions reduced to : ', num2str(length(reduced_solutions.array)) ' for misorietation of CPPs  < ',num2str(theta_CPP_max),'°'] );
            case 4
                reduced_solutions = Solution_array( Slip_solution(), reduced_solutions, handles.cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h');
                updateLog_MartCalc(hObject, handles, ['Solutions reduced to : ', num2str(length(reduced_solutions.array)) ' for habit plane misorientation to CP-planes  < ',num2str(theta_n_max),'°'] );
            case 5
                reduced_solutions = Solution_array( Slip_solution(), reduced_solutions, 'det', delta_determinant_max,  det(martensite.U));
                updateLog_MartCalc(hObject, handles, ['Solutions reduced to : ', num2str(length(reduced_solutions.array)) ' for (non-physical) volume change  < ',num2str(delta_determinant_max)] );
            case 6
                reduced_solutions = Solution_array( Slip_solution(), reduced_solutions, handles.KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );
                updateLog_MartCalc(hObject, handles, ['Solutions reduced to : ', num2str(length(reduced_solutions.array)) ' for a maximum deviation angle of KS-directions  < ',num2str(theta_KS_max),'°'] );
            case 7
                reduced_solutions = Solution_array( Slip_solution(), reduced_solutions, handles.NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);
                updateLog_MartCalc(hObject, handles, ['Solutions reduced to : ', num2str(length(reduced_solutions.array)) ' for a maximum deviation angle of NW-directions  < ',num2str(theta_NW_max),'°'] );
        end
    end
    handles.reduced_solutions = reduced_solutions;
    updateLog_MartCalc(hObject, handles, 'Filtering of IPS solutions after specified criteria completed.');
    guidata(hObject, handles);
else
    updateLog_MartCalc(hObject, handles, 'Note: No filters for solutions specified.');
end
