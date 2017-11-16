function [ ] = delete_asc_MartCalc( hObject, eventdata ) % here both: hObject and eventdata are necessary that it works!!!
% DELETE_ASC_MARTCALC 

% PET 15.11.17 - handles.pan_asc ---> handles.(handles.paneltag) ---> main_asc_panel

% retrieve data of GUI
handles = guidata(hObject);
button_parent = get(hObject,'Parent'); % get parent panel of asc
main_asc_panel = get(button_parent,'Parent'); % get main panel for asc

% children(3) is the textfield with a description of the actual criterion
% its .Value property is the "internal number" for our sorting of selection
% criteria
criterion_number = button_parent.Children(3).Value; % children(3) is the textfield with a description of the actual criterion

% PET 16.11.17
switch main_asc_panel
    case handles.pan_asc_IPS
        status     = handles.asc_status_IPS;
        asc_number = handles.asc_number_IPS;
        asc_list = handles.asc_list_IPS;
    case handles.pan_asc_ILS
        % criteria for ILS solutions
        status     = handles.asc_status_ILS;
        asc_number = handles.asc_number_ILS;
        asc_list = handles.asc_list_ILS;
    otherwise
        error('This parent is not considered for deletion of panels')
end


%% 
postion_in_asc_list = status(criterion_number); %handles.asc_status(criterion_number); 
% delete actual asc panel
% NOTE: with deletion of this panel, the following panels in the list are
% switched one position forward in the list of .Children!
% % % handles.pan_asc.Children(handles.asc_number-postion_in_asc_list+1).delete(); 
main_asc_panel.Children(size(main_asc_panel.Children,1)+1-postion_in_asc_list).delete() 
asc_number = asc_number - 1; % decrease number of asc
status(criterion_number) = 0; % set criterion inactive


% update positions of the following list entries
% % while(handles.asc_list(postion_in_asc_list) > 0)
while(postion_in_asc_list <= size(main_asc_panel.Children,1))
   % update position of the panel
      main_asc_panel.Children(size(main_asc_panel.Children,1)+1-postion_in_asc_list).Position(2)=...
       main_asc_panel.Children(size(main_asc_panel.Children,1)+1-postion_in_asc_list).Position(2)+ ...
       ((main_asc_panel.Position(4)/8)/main_asc_panel.Position(4));
   
   % update number of position in list
   main_asc_panel.Children(size(main_asc_panel.Children,1)+1-postion_in_asc_list).Children(4).String = num2str(postion_in_asc_list);
   
   % update status of criterion in order to fit new position in list
   status(main_asc_panel.Children(size(main_asc_panel.Children,1)+1-postion_in_asc_list).Children(3).Value) = postion_in_asc_list; 
   
   asc_list(postion_in_asc_list) = main_asc_panel.Children(postion_in_asc_list).Children(3).Value;
   
   postion_in_asc_list = postion_in_asc_list+1;
end
%%


switch main_asc_panel
    case handles.pan_asc_IPS
        handles.asc_status_IPS = status;
        handles.asc_number_IPS = asc_number;
        handles.asc_list_IPS = asc_list;
    case handles.pan_asc_ILS
        % criteria for ILS solutions
        handles.asc_status_ILS = status;
        handles.asc_number_ILS = asc_number;
        handles.asc_list_ILS = asc_list;
end

% update GUI data
guidata(main_asc_panel, handles); % same as - 
%guidata(hobject, handles);
end

