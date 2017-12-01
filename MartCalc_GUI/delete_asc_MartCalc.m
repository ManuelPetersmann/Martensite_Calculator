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

%criterion_number = button_parent.Children(3).Value; % children(3) is the textfield with a description of the actual criterion
butt_par_child = get(button_parent,'Children'); 
criterion_number = get(butt_par_child(3),'Value');

%criterion_number
%class('criterion_number')
%criterion_number = str2num(criterion_number)

% PET 16.11.17 - Note: At the end of this file the variables are written back =)
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
    case handles.pan_asc_blocks
        status     = handles.asc_status_blocks;
        asc_number = handles.asc_number_blocks;
        asc_list = handles.asc_list_blocks;
    otherwise
        error('This parent is not considered for deletion of panels')
end


%% 
postion_in_asc_list = status(criterion_number); %handles.asc_status(criterion_number); 
% delete actual asc panel
% NOTE: with deletion of this panel, the following panels in the list are
% switched one position forward in the list of .Children!

main_asc_pan_childs = get(main_asc_panel,'Children');% PET 30.11.17

% % % handles.pan_asc.Children(handles.asc_number-postion_in_asc_list+1).delete(); 
delete( main_asc_pan_childs( size(main_asc_pan_childs,1)+1-postion_in_asc_list ) ) %.delete() 
asc_number = asc_number - 1; % decrease number of asc
status(criterion_number) = 0; % set criterion inactive

main_asc_pan_childs = get(main_asc_panel,'Children');% PET 30.11.17
main_asc_pan_pos =    get(main_asc_panel,'Position');

% size(main_asc_pan_childs,1)
% size(main_asc_panel.Children,1)
% postion_in_asc_list

% update positions of the following list entries
% % while(handles.asc_list(postion_in_asc_list) > 0)
while(postion_in_asc_list <= size(main_asc_pan_childs,1))
    
   %% example i got from askign mathworks 
   %pos = get(handles.panel.children(4), 'Position');
   %pos(3) = some value;
   %set(handles.panel.Children(4), 'Position', pos);
   %% 
   %disp('lol')
    
   % update position of the panel
   main_asc_pan_childs = get(main_asc_panel,'Children');
   index = size(main_asc_pan_childs,1) +1 - postion_in_asc_list; 
   %temp_position = get( main_asc_pan_childs( index ) , 'Position');
   temp_position = get( main_asc_pan_childs( index ) , 'Position');
   temp_position(2) = temp_position(2) + (main_asc_pan_pos(4)/8)/main_asc_pan_pos(4);
   %set( temp_position(2), ...    main_asc_pan_childs( index ).Position(2) = ...
   %    temp_position(2) +  ( main_asc_pan_pos(4)/8) / main_asc_pan_pos(4) );
   %   main_asc_pan_childs( size(main_asc_pan_childs,1)+1-postion_in_asc_list ).Position(2)+ ...
   %   ((main_asc_pan_pos(4)/8)/main_asc_pan_pos(4));
   set(main_asc_pan_childs( index ), 'Position', temp_position);  

   
   child_childs = get(main_asc_pan_childs(index),'Children');
   % update number of position in list
   %main_asc_pan_childs(index).Children(4).String = num2str(postion_in_asc_list);
   set(child_childs(4),'String',num2str(postion_in_asc_list) );
   
   % update status of criterion in order to fit new position in list
   status( get(child_childs(3),'Value') ) = postion_in_asc_list; 
   %set( status( get( child_childs(3),'Value' ) ) , postion_in_asc_list); 
   
   bla = get(main_asc_pan_childs(postion_in_asc_list),'Children');
   %asc_list(postion_in_asc_list) = main_asc_pan_childs(postion_in_asc_list).Children(3).Value;
   asc_list(postion_in_asc_list) = get( bla(3),'Value');
   
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
    case handles.pan_asc_blocks
        handles.asc_status_blocks = status;
        handles.asc_number_blocks = asc_number;
        handles.asc_list_blocks = asc_list;
end

% update GUI data
guidata(main_asc_panel, handles); % same as - 
%guidata(hobject, handles);
end

