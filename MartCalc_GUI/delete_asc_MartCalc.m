function [ ] = delete_asc_MartCalc( hObject,eventdata)
%DELETE_ASC_MARTCALC Summary of this function goes here
%   Detailed explanation goes here

% retrieve data of GUI
handles = guidata(hObject);
button_parent = get(hObject,'Parent'); % get parent panel of asc

main_asc_panel = get(button_parent,'Parent'); % get main panel for asc

% children(3) is the textfield with a description of the actual criterion
% its .Value property is the "internal number" for our sorting of selection
% criteria
criterion_number = button_parent.Children(3).Value; % children(3) is the textfield with a description of the actual criterion

% EHL: add functionality for updating the asc list in the GUI!!!
% ....
postion_in_asc_list = handles.asc_status(criterion_number);
% delete actual asc panel
% NOTE: with deletion of this panel, the following panels in the list are
% switched one position forward in the list of .Children!
% % % handles.pan_asc.Children(handles.asc_number-postion_in_asc_list+1).delete(); 
handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-postion_in_asc_list).delete(); 
handles.asc_number = handles.asc_number - 1; % decrease number of asc
handles.asc_status(criterion_number) = 0; % set criterion inactive


% update positions of the following list entries
% % while(handles.asc_list(postion_in_asc_list) > 0)
while(postion_in_asc_list <= size(handles.pan_asc.Children,1))
   % update position of the panel
   handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-postion_in_asc_list).Position(2)=...
       handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-postion_in_asc_list).Position(2)+((handles.pan_asc.Position(4)/8)/handles.pan_asc.Position(4));
   
   % update number of position in list
   handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-postion_in_asc_list).Children(4).String = num2str(postion_in_asc_list);
   
   % update status of criterion in order to fit new position in list
   handles.asc_status(handles.pan_asc.Children(size(handles.pan_asc.Children,1)+1-postion_in_asc_list).Children(3).Value) = postion_in_asc_list; 
   
   handles.asc_list(postion_in_asc_list)=handles.pan_asc.Children(postion_in_asc_list).Children(3).Value;
   
   postion_in_asc_list = postion_in_asc_list+1;
end

% update GUI data
% % % guidata(Martensite_Calculator, handles);
guidata(main_asc_panel, handles);
end

