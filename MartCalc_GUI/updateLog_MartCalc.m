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

log = get(handles.log_lb, 'string'); % fetch log
nol_log = size(log,1); % Number of Lines - in actual log (size of log)

% first try: of wrapping text automatically in log:
% if length(str_log_update) > 67
%     str_log_update = [str_log_update(1:67) sprintf('\n') str_log_update(68:length(str_log_update))];
% end
% log{nol_log+1,1} = [time,' - ',str_log_update]; % extend log
% set(handles.log_lb, 'string', log, 'value', 1); % update log
% guidata(hObject, handles); % Update handles structure
%end
%
% second try: wrapped text with function
% handles.log_lb - uicontrols object - ok!
[outstring, ~] = textwrap(handles.log_lb, {[time,' - ',str_log_update] } ); % create new string for log


% old extend log
%log{nol_log+1,1} = outstring %[time,' - ',str_log_update]; 
%
% new extend log
log = [log; outstring]; % concatenating cell arrays

% update log
set(handles.log_lb, 'string', log, 'value', nol_log +1 ); 
%currView = get(handles.log_lb,'ListBoxTop'); 
%set(handles.log_lb, 'Listboxtop',topIndex);
%set(handles.log_lb, 'Listboxtop', currView); 

guidata(handles.log_lb, handles); % Update handles structure

drawnow(); % update GUI
end

