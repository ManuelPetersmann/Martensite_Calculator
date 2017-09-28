function varargout = Martensite_Calculator(varargin)
%MARTENSITE_CALCULATOR MATLAB code file for Martensite_Calculator.fig
%      MARTENSITE_CALCULATOR, by itself, creates a new MARTENSITE_CALCULATOR or raises the existing
%      singleton*.
%
%      H = MARTENSITE_CALCULATOR returns the handle to a new MARTENSITE_CALCULATOR or the handle to
%      the existing singleton*.
%
%      MARTENSITE_CALCULATOR('Property','Value',...) creates a new MARTENSITE_CALCULATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Martensite_Calculator_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MARTENSITE_CALCULATOR('CALLBACK') and MARTENSITE_CALCULATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MARTENSITE_CALCULATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Martensite_Calculator

% Last Modified by GUIDE v2.5 13-Sep-2017 16:38:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Martensite_Calculator_OpeningFcn, ...
                   'gui_OutputFcn',  @Martensite_Calculator_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Martensite_Calculator is made visible.
function Martensite_Calculator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Martensite_Calculator
handles.output = hObject;

% initialize array for keeping track of the active selection criteria
% 28.08.2017: currently there are 7 possible criteria for selection of solutions
% default at start of MartCalc-GUI: all = 0 --> inactive
% asc = Active Selection Criteria
handles.asc_status = zeros(1,7);
handles.asc_number = 0;
handles.asc_list = zeros(1,7);

handles.log_status = 0; % variable for check if log has already been changed for a first time

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes Martensite_Calculator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Martensite_Calculator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_MartCalc.
function start_MartCalc_Callback(hObject, eventdata, handles)
% hObject    handle to start_MartCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clear log and start from scratch
% % % start_time = datestr(now,'HH:MM:SS.FFF'); % get time
% % % new_log = {[start_time,' - ','Procedure started.']}; % start new log as cell array for several lines
% % % set(handles.log_lb, 'string', new_log, 'value', 1); % update log
updateLog_MartCalc(hObject, handles, 'Procedure started.')

%% Set up coordinate systems and lattice parameters
martensite = Martensite(); % creates martensite object
austenite = Base();

% read user input from GUI for determination of solutions
updateLog_MartCalc(hObject, handles, 'Retrieving input from GUI')
get_input_MartCalc;

% calculate bain strains and correspondence matrix based on given
% lattice constants
% % % updateLog_MartCalc(hObject, handles, 'Calculating Bain strain and correspondence matrix.')
% % % Bain_and_Correspondence_mavalX12;

% start determination of possible solution and filtering according to
% active selection criteria
calc_mech = handles.popup_calc_mech.Value;

updateLog_MartCalc(hObject, handles, 'Determination of solutions started.')
switch calc_mech
    case 1
        updateLog_MartCalc(hObject, handles, 'maraging_block_sym_doubleshear - run')
        maraging_block_sym_doubleshear; 
    case 2
        updateLog_MartCalc(hObject, handles, 'maraging_block_sym_doubleshear_specific_slipsys_test - run')
        maraging_block_sym_doubleshear_specific_slipsys_test;
    case 3
        updateLog_MartCalc(hObject, handles, 'maraging_MarescaCurtin_test - run')
        maraging_MarescaCurtin_test;
    case 4
        updateLog_MartCalc(hObject, handles, 'maraging_multiple_shears - run')
        maraging_multiple_shears;
    case 5
        updateLog_MartCalc(hObject, handles, 'maraging_variable_doubleshear - run')
        maraging_variable_doubleshear;
end
updateLog_MartCalc(hObject, handles, 'Determination and filtering of solutions completed.')

% --- Executes on selection change in lsc_popup.
function lsc_popup_Callback(hObject, eventdata, handles)
% hObject    handle to lsc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lsc_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lsc_popup

switch hObject.Value
    case 1 % Criterion 1 has been chosen: Minimum slip plane density
        if handles.asc_status(1) == 0
            handles.asc_number = handles.asc_number + 1; % increase number of asc
            
            criterion_name = 'Minimum slip plane density';
            default_value = 10.0;
            
            handles.asc_list(handles.asc_number) = hObject.Value; % keep track of which criterion is at which point in the asc list
            handles.asc_status(1) = handles.asc_number; % set = number in row, in order to show that crit is already active and when it is to be applied
            
            handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, hObject.Value);
            
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, 'Criterion - "Minimum slip plane density" is already active!.')
        end
    %
    case 2 % Criterion 2 has been chosen: Maximum shape strain
        if handles.asc_status(2) == 0
            handles.asc_number = handles.asc_number + 1; % increase number of asc
            
            criterion_name = 'Maximum shape strain';
            default_value = 0.4;
            
            handles.asc_list(handles.asc_number) = hObject.Value; % keep track of which criterion is at which point in the asc list
            handles.asc_status(2) = handles.asc_number; % set = number in row, in order to show that crit is already active and when it is to be applied
            
            handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, hObject.Value);
           
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, 'Criterion - "Maximum shape strain" is already active!.')
        end
    case 3 % Criterion 3 has been chosen: Maximum misorientation of CPPs {110}_alpha and {111}_gamma
        if handles.asc_status(3) == 0
            handles.asc_number = handles.asc_number + 1; % increase number of asc
            
            criterion_name = 'Maximum misorientation of CPPs {110}_alpha and {111}_gamma';
            default_value = 3.0;
            
            handles.asc_list(handles.asc_number) = hObject.Value; % keep track of which criterion is at which point in the asc list
            handles.asc_status(3) = handles.asc_number; % set = number in row, in order to show that crit is already active and when it is to be applied
            
            handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, hObject.Value);
                                      
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, 'Criterion - "Maximum misorientation of CPPs {110}_alpha and {111}_gamma" is already active!.')
        end
    case 4 % Criterion 4 has been chosen: Maximum misorientation of block HP to {111}_gamma
        if handles.asc_status(4) == 0
            handles.asc_number = handles.asc_number + 1; % increase number of asc
            
            criterion_name = 'Maximum misorientation of block HP to {111}_gamma';
            default_value = 5.0;
            
            handles.asc_list(handles.asc_number) = hObject.Value; % keep track of which criterion is at which point in the asc list
            handles.asc_status(4) = handles.asc_number; % set = number in row, in order to show that crit is already active and when it is to be applied
            
            handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, hObject.Value);
                                      
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, 'Criterion - "Maximum misorientation of block HP to {111}_gamma" is already active!.')
        end
    case 5 % Criterion 5 has been chosen: Maximum deviation of determinant det(F) of transformation
        if handles.asc_status(5) == 0
            handles.asc_number = handles.asc_number + 1; % increase number of asc
            
            criterion_name = 'Maximum deviation of determinant det(F) of transformation';
            default_value = 0.0001;
            
            handles.asc_list(handles.asc_number) = hObject.Value; % keep track of which criterion is at which point in the asc list
            handles.asc_status(5) = handles.asc_number; % set = number in row, in order to show that crit is already active and when it is to be applied
            
            handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, hObject.Value);
            
            guidata(hObject, handles); % Update handles structure
        else
        	updateLog_MartCalc(hObject, handles, 'Criterion - "Maximum deviation of determinant det(F) of transformation" is already active!.')
        end
    case 6 % Criterion 6 has been chosen: Maximum deviation from KS OR
        if handles.asc_status(6) == 0
            handles.asc_number = handles.asc_number + 1; % increase number of asc
            
            criterion_name = 'Maximum deviation from KS OR';
            default_value = 10.0;
            
            handles.asc_list(handles.asc_number) = hObject.Value; % keep track of which criterion is at which point in the asc list
            handles.asc_status(6) = handles.asc_number; % set = number in row, in order to show that crit is already active and when it is to be applied
            
            handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, hObject.Value);
            
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, 'Criterion - "Maximum deviation from KS OR" is already active!.')
        end
    case 7 % Criterion 7 has been chosen: Maximum deviation from NW OR
        if handles.asc_status(7) == 0
            handles.asc_number = handles.asc_number + 1; % increase number of asc
            
            criterion_name = 'Maximum deviation from NW OR';
            default_value = 10.0;
            
            handles.asc_list(handles.asc_number) = hObject.Value; % keep track of which criterion is at which point in the asc list
            handles.asc_status(7) = handles.asc_number; % set = number in row, in order to show that crit is already active and when it is to be applied
            
            handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, hObject.Value);
            
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, 'Criterion - "Maximum deviation from NW OR" is already active!.')
        end
end
% % % guidata(hObject, handles); % Update handles structure


% --- Executes during object creation, after setting all properties.
function lsc_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lsc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssf_edtxt_sp_val_Callback(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sp_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssf_edtxt_sp_val as text
%        str2double(get(hObject,'String')) returns contents of ssf_edtxt_sp_val as a double


% --- Executes during object creation, after setting all properties.
function ssf_edtxt_sp_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sp_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssf_edtxt_sd_val_Callback(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sd_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssf_edtxt_sd_val as text
%        str2double(get(hObject,'String')) returns contents of ssf_edtxt_sd_val as a double


% --- Executes during object creation, after setting all properties.
function ssf_edtxt_sd_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sd_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssf_edtxt_sp_name_Callback(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sp_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssf_edtxt_sp_name as text
%        str2double(get(hObject,'String')) returns contents of ssf_edtxt_sp_name as a double


% --- Executes during object creation, after setting all properties.
function ssf_edtxt_sp_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sp_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssf_edtxt_sd_name_Callback(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sd_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssf_edtxt_sd_name as text
%        str2double(get(hObject,'String')) returns contents of ssf_edtxt_sd_name as a double


% --- Executes during object creation, after setting all properties.
function ssf_edtxt_sd_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssf_edtxt_sd_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_edtxt_aust_val_Callback(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_aust_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_edtxt_aust_val as text
%        str2double(get(hObject,'String')) returns contents of lc_edtxt_aust_val as a double


% --- Executes during object creation, after setting all properties.
function lc_edtxt_aust_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_aust_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_edtxt_mart_val_Callback(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_mart_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_edtxt_mart_val as text
%        str2double(get(hObject,'String')) returns contents of lc_edtxt_mart_val as a double


% --- Executes during object creation, after setting all properties.
function lc_edtxt_mart_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_mart_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_edtxt_aust_name_Callback(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_aust_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_edtxt_aust_name as text
%        str2double(get(hObject,'String')) returns contents of lc_edtxt_aust_name as a double


% --- Executes during object creation, after setting all properties.
function lc_edtxt_aust_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_aust_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_edtxt_mart_name_Callback(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_mart_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_edtxt_mart_name as text
%        str2double(get(hObject,'String')) returns contents of lc_edtxt_mart_name as a double


% --- Executes during object creation, after setting all properties.
function lc_edtxt_mart_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_edtxt_mart_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_main_Callback(hObject, eventdata, handles)
% hObject    handle to menu_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in log_lb.
function log_lb_Callback(hObject, eventdata, handles)
% hObject    handle to log_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns log_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from log_lb


% --- Executes during object creation, after setting all properties.
function log_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to log_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent11_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent11 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent11 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent12_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent12 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent12 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent13_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent13 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent13 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent21_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent21 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent21 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent22_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent22 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent22 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent23_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent23 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent23 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent31_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent31 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent31 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent32_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent32 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent32 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_parent33_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_parent33 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_parent33 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_parent33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_parent33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product11_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product11 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product11 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product12_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product12 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product12 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product13_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product13 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product13 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product21_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product21 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product21 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product22_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product22 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product22 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product23_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product23 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product23 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product31_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product31 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product31 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product32_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product32 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product32 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basevec_edtxt_product33_Callback(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basevec_edtxt_product33 as text
%        str2double(get(hObject,'String')) returns contents of basevec_edtxt_product33 as a double


% --- Executes during object creation, after setting all properties.
function basevec_edtxt_product33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basevec_edtxt_product33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt11_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt11 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt11 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt12_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt12 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt12 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt13_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt13 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt13 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt21_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt21 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt21 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt22_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt22 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt22 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt23_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt23 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt23 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt31_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt31 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt31 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt32_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt32 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt32 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corrmat_edtxt33_Callback(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrmat_edtxt33 as text
%        str2double(get(hObject,'String')) returns contents of corrmat_edtxt33 as a double


% --- Executes during object creation, after setting all properties.
function corrmat_edtxt33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrmat_edtxt33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_calc_mech.
function popup_calc_mech_Callback(hObject, eventdata, handles)
% hObject    handle to popup_calc_mech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_calc_mech contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_calc_mech


% --- Executes during object creation, after setting all properties.
function popup_calc_mech_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_calc_mech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
