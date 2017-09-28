function [] = updateLog_MartCalc(hObject, handles, str_log_update)
%UPDATE_GUI_LOG Summary of this function goes here
%   Detailed explanation goes here
if handles.log_status == 0
    time = datestr(now,'HH:MM:SS.FFF'); % get time
    log = {[time,' - ',str_log_update]}; % start new log
	set(handles.log_lb, 'string', log, 'value', 1); % update log
    
    handles.log_status = 1;
    
    guidata(hObject, handles); % Update handles structure
else
    log = get(handles.log_lb, 'string'); % fetch log 
    nol_log = size(log,1); % number of lines in actual log
    time = datestr(now,'HH:MM:SS.FFF'); % get time
    log{nol_log+1,1} = [time,' - ',str_log_update]; % extend log
    set(handles.log_lb, 'string', log, 'value', 1); % update log
    
    guidata(hObject, handles); % Update handles structure
end

drawnow(); % update GUI
end

