%% selectin criteria
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
    theta_p_max = str2num(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-handles.asc_status(3)).Children(2).String);
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