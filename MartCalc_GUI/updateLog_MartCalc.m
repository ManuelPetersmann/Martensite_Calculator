function [] = updateLog_MartCalc(hObject, handles, str_log_update)  
%call updateLog_MartCalc(handles, str_log_update)

time = datestr(now,'HH:MM:SS.FFF'); % get time

% PET: is this necessary?
% if handles.log_status == 0
%     log = {[time,' - ',str_log_update]}; % start new log = cell-array
%     set(handles.log_lb, 'string', log, 'value', 1); % update log
%     handles.log_status = 1;
%     guidata(hObject, handles); % Update handles structure
% else

log = get(handles.log_lb, 'string') % fetch log
nol_log = size(log,1); % Number of Lines - in actual log (size of log)

% first try: of wrapping text automatically in log:
%handles.log_lb.Position(3) % width
% handles.log_lb.Position(4) % height
%if length(str_log_update) > 67
%    str_log_update = [str_log_update(1:67) newline str_log_update(68:length(str_log_update))];
%end

% TODO TODO TODO 

% second try: wrapped text with function
[outstring, newpos] = textwrap(handles.log_lb, {time,' - ',str_log_update } ) % create new string for log

%log{nol_log+1,1} = outstring %[time,' - ',str_log_update]; % extend log
log = [log; outstring]

set(handles.log_lb, 'String',log); %,'Position',newpos)

% old version
%set(handles.log_lb, 'string', log, 'value', 1); % update log

% set log to last message
%set(handles.listbox1, 'Listboxtop', nol_log); % set log to last position

guidata(handles.log_lb, handles); % Update handles structure
%end

drawnow(); % update GUI
end

