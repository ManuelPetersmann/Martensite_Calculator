function [] = updateLog_MartCalc(hObject, handles, str_log_update)
%UPDATE_GUI_LOG Summary of this function goes here
%   Detailed explanation goes here

time = datestr(now,'HH:MM:SS.FFF'); % get time

% if handles.log_status == 0
%     log = {[time,' - ',str_log_update]}; % start new log = cell-array
%     set(handles.log_lb, 'string', log, 'value', 1); % update log
%     handles.log_status = 1;
%     guidata(hObject, handles); % Update handles structure
% else
log = get(handles.log_lb, 'string'); % fetch log
nol_log = size(log,1); % number of messages in actual log
if length(str_log_update) > 67
    str_log_update = [str_log_update(1:67) newline str_log_update(68:length(str_log_update))];
end
log{nol_log+1,1} = [time,' - ',str_log_update]; % extend log
set(handles.log_lb, 'string', log, 'value', 1); % update log
guidata(hObject, handles); % Update handles structure
%end

drawnow(); % update GUI
end

