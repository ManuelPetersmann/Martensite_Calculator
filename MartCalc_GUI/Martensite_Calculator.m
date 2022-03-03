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

% Last Modified by GUIDE v2.5 22-Sep-2021 15:53:58

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

% Initialise tabs
% PET 23.10.17: Commented for one page GUI for now - uncomment to enable tabbing
% handles.tabManager = TabManager( hObject );
% % Set-up a selection changed function on the create tab groups
% tabGroups = handles.tabManager.TabGroups;
% for tgi=1:length(tabGroups)
%     set(tabGroups(tgi),'SelectionChangedFcn',@tabChangedCB)
% end

% initialize array for keeping track of the active selection criteria
% 28.08.2017: currently there are 7 possible criteria for selection of IPS solutions
% default at start of MartCalc-GUI: all = 0 --> inactive
% asc = Active Selection Criteria
handles.asc_status_IPS = zeros(1,7);
handles.asc_number_IPS = 0;
handles.asc_list_IPS = zeros(1,7);
handles.asc_IPS_maxcrit = 7;
%
% criteria for ILS solutions
handles.asc_status_ILS = zeros(1,7);
handles.asc_number_ILS = 0;
handles.asc_list_ILS = zeros(1,7);
handles.asc_ILS_maxcrit = 7;

handles.asc_status_blocks = zeros(1,3);
handles.asc_number_blocks = 0;
handles.asc_list_blocks = zeros(1,3);
handles.asc_blocks_maxcrit = 3;

% create austenite and martensite objects
handles.martensite = Martensite(); % creates martensite object
handles.austenite  = Base();

% integrated into Bravais object: CPPs-close packed planes and CP_dirs - close packed directions 
handles.austenite.CPPs    = all_from_family_perms( [1 1 1] ); % close packed planes of gamma-lattice - formerly 'cpps_gamma'
handles.austenite.CP_dirs = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0] - formerly 'KS'
handles.NW = all_from_family_perms( [1 2 1], false );

handles.input_status   = true; % will be set to false if something is wrong with the input
%handles.log_status = 0; % variable for check if log has already been changed for a first time

handles.InterfaceObj = findobj(handles.figure1,'Style','pushbutton','-or','Style','popupmenu'); % variable to disable interface during calculations
%handles.InterfaceObj = findobj(handles.figure1,'Style','pushbutton','-or','Tag','popup_sorting_ILS','-or','Tag','popup_sorting_IPS'); % variable to disable interface during calculations
%handles.InterfaceObj = findobj(handles.figure1,'Enable','on','-and','-not','Tag','filename_results_edittext'); % variable to disable interface during calculations

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Martensite_Calculator_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;



%%#########################################################################
%% LATH PART %%

% --- Executes on selection change in popup_calc_lath_level.
function popup_calc_lath_level_Callback(hObject, eventdata, handles)
% hObject    handle to popup_calc_lath_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popup_calc_lath_level contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_calc_lath_level

lath_type = get(hObject,'Value');
switch lath_type %handles.popup_calc_lath_level.Value
    case 1 % or or or...
        set(handles.pan_lsc_IPS,'visible','on')
        set(handles.pan_asc_IPS,'visible','on')
        set(handles.pan_sort_IPS,'visible','on')
        set(handles.pan_lsc_ILS,'visible','off')
        set(handles.pan_asc_ILS,'visible','off')
        set(handles.pan_sort_ILS,'visible','off')
    case 2 % or or or
        set(handles.pan_lsc_ILS,'visible','on')
        set(handles.pan_asc_ILS,'visible','on')
        set(handles.pan_sort_ILS,'visible','on')
        set(handles.pan_lsc_IPS,'visible','off')
        set(handles.pan_asc_IPS,'visible','off')
        set(handles.pan_sort_IPS,'visible','off')
end
guidata(hObject, handles);


%% --- Executes on button press in pushbutton_start_lath_calc.
function pushbutton_start_lath_calc_Callback(hObject, eventdata, handles)

% disable interface during function call
set(handles.InterfaceObj,'Enable','off');
% enable interface again
%set(handles.InterfaceObj,'Enable','on');

updateLog_MartCalc(hObject, handles, '---------------------------------------------')
updateLog_MartCalc(hObject, handles, 'Retrieving input from GUI')
% read user input from GUI for determination of solutions
get_input_MartCalc;

if handles.input_status
    % start determination of possible solution and subsequent filtering according to active selection criteria
    % store solution objects in an object array
    switch handles.popup_calc_lath_level.Value
        case 1
            calculation_method = 'IPS by doubleshear shear, incremental minimization of IPS condition abs(lambda2 -1.) at the lath level';
            updateLog_MartCalc(hObject, handles, [calculation_method,' - started, calculating...']);
            handles.martensite.IPS_solutions.calculation_method = calculation_method;
            handles.martensite.IPS_solutions = doubleshear_variable_shear_mags(handles.martensite, handles.austenite);
            %
            updateLog_MartCalc(hObject, handles, ['Determination of IPS solutions for laths completed: ' num2str(size(handles.martensite.IPS_solutions.array,2)),' solutions found.'] );
            % ILS
        case 2
            calculation_method = 'ILS, doubleshear, incremental minimization of unrotated vector residual';
            updateLog_MartCalc(hObject, handles, [calculation_method,' - started, calculating...']);
            handles.martensite.ILS_solutions.calculation_method = calculation_method;
            %handles.martensite.ILS_solutions = invariant_line_strain_fixed_increment(handles.martensite, handles.austenite); 
            handles.martensite.ILS_solutions = invariant_line_strain_line_search(handles.martensite, handles.austenite); 
            %
            updateLog_MartCalc(hObject, handles, ['Determination of ILS solutions for laths completed: ' num2str(size(handles.martensite.ILS_solutions.array,2)),' solutions found.'] );
            
            %% other cases can be added here         
            %         updateLog_MartCalc(hObject, handles, 'multiple shears incremental minimization - started')
            %         maraging_multiple_shears;
            %     case 3
            %         updateLog_MartCalc(hObject, handles, '_MarescaCurtin_test - run')
            %         maraging_MarescaCurtin_test;
    end
    
    %% filter solutions for volume change not necessary - volume does not change in my calculations!
%     l1 = length(handles.martensite.IPS_solutions.array);
%     %  formerly - case 8 - now default reduction to this value!
%     delta_determinant_max = 0.001; % maximum 0.1% non-physical volume change
%     updateLog_MartCalc(hObject, handles, ['checking solutions for non-physical volume change > ',num2str(delta_determinant_max*100),'%...'] );
%     % handles.reduced_solutions = 
%     handles.martensite.IPS_solutions = Solution_array( IPS_solution(), handles.martensite.IPS_solutions, 'delta_determinant_max', delta_determinant_max,  det(handles.martensite.U) );
%     l2 = length(handles.martensite.IPS_solutions.array);
%     if l1 ~= l2
%         updateLog_MartCalc(hObject, handles, [num2str(l2-l1),' solutions discarded' ] );
%     else
%         updateLog_MartCalc(hObject, handles,'all solutions are valid');
%     end
    update_lath_selection_criteria;
else
    updateLog_MartCalc(hObject, handles, 'Calculation could not be started - insufficient input - see above log messages.');
end
guidata(hObject, handles);

% enable interface again
set(handles.InterfaceObj,'Enable','on');



%%#########################################################################
%% --- Executes on selection change in lsc_popup_IPS.
function lsc_popup_IPS_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns lsc_popup_IPS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lsc_popup_IPS

selected_crit = get(hObject,'Value');
switch selected_crit
    case 1 % Criterion 1 has been chosen: Minimum slip plane density
        criterion_name = 'Minimum average transformation-dislocation spacing (stepwidth on interface).';
        if handles.asc_status_IPS(1) == 0 % if inactive
            handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc
            default_value = 5.0;
            handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_IPS(1) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_IPS, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
    %
    case 2 % Criterion 2 has been chosen: Maximum shape strain
        criterion_name = 'Maximum (total) shape strain (eps_ips) of invariant plane strain.';
        if handles.asc_status_IPS(2) == 0
            handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc
            default_value = 0.6;
            handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_IPS(2) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_IPS, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
    case 3 % Criterion 3 has been chosen: Maximum misorientation of CPPs {110}_alpha and {111}_gamma
        criterion_name = 'Maximum misorientation of {111}_gamma to {110}_alpha';
        if handles.asc_status_IPS(3) == 0
            handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc
            default_value = 2.0;
            handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_IPS(3) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_IPS, criterion_name, default_value, selected_crit );                        
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
    case 4 % Criterion 4 has been chosen: Maximum misorientation of block HP to {111}_gamma
        criterion_name = 'Maximum misorientation of invariant plane to {111}_gamma.';
        if handles.asc_status_IPS(4) == 0
            handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc
            default_value = 20.0;
            handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_IPS(4) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_IPS, criterion_name, default_value, selected_crit );                       
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
    case 5 % Criterion 6 has been chosen: Maximum deviation from KS OR
        criterion_name = 'Maximum deviation of perfect KS direction parallelism.';
        if handles.asc_status_IPS(5) == 0
            handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc
            default_value = 5.0;
            handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_IPS(5) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_IPS, criterion_name, default_value, selected_crit);
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
    case 6 % Criterion 7 has been chosen: Maximum deviation from NW OR
        criterion_name = 'Maximum deviation from perfect NW direction parallelism';
        if handles.asc_status_IPS(6) == 0
            handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc    
            default_value = 8.0;
            handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_IPS(6) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied       
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_IPS, criterion_name, default_value, selected_crit); 
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
    case 7 % Criterion 8 has been chosen: Maximum deviation of preferred invariant line to invariant habit plane
        criterion_name = 'Maximum tolerance angle between preferred invariant line and habit plane';
        if handles.asc_status_IPS(7) == 0
            handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc
            default_value = 3.0;
            handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_IPS(7) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_IPS, criterion_name, default_value, selected_crit);
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
%     case 8 % Criterion 5 has been chosen: Maximum deviation of determinant det(F) of transformation
%         if handles.asc_status_IPS(8) == 0
%             handles.asc_number_IPS = handles.asc_number_IPS + 1; % increase number of asc
%             criterion_name = 'Maximum deviation of theoretical volume change from Bain strain.';
%             default_value = 0.001;
%             handles.asc_list_IPS(handles.asc_number_IPS) = selected_crit; % keep track of which criterion is at which point in the asc list
%             handles.asc_status_IPS(8) = handles.asc_number_IPS; % set = number in row, in order to show that crit is already active and when it is to be applied
%             handles = create_asc_panel_MartCalc(handles, criterion_name, default_value, selected_crit);
%             guidata(hObject, handles); % Update handles structure
%         else
%             updateLog_MartCalc(hObject, handles, 'Criterion - "Maximum deviation of theoretical volume change from Bain strain" is already active!')
%         end
end


% --- Executes on selection change in lsc_popup_ILS.
function lsc_popup_ILS_Callback(hObject, eventdata, handles)

selected_crit = get(hObject,'Value');
switch selected_crit
    case 1 % Criterion 1 : Minimum slip plane density
        criterion_name = 'Minimum average transformation-dislocation spacing (stepwidth on interface).';
        if handles.asc_status_ILS(1) == 0 % if inactive
            handles.asc_number_ILS = handles.asc_number_ILS + 1; % increase number of asc
            default_value = 5.0;
            handles.asc_list_ILS(handles.asc_number_ILS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_ILS(1) = handles.asc_number_ILS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_ILS, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
        %
    case 2 % Criterion 2 : Maximum misorientation of CPPs {110}_alpha and {111}_gamma
        criterion_name = 'Maximum misorientation of {111}_gamma to {110}_alpha';
        if handles.asc_status_ILS(2) == 0
            handles.asc_number_ILS = handles.asc_number_ILS + 1; % increase number of asc
            default_value = 2.0;
            handles.asc_list_ILS(handles.asc_number_ILS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_ILS(2) = handles.asc_number_ILS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_ILS, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
        %
    case 3 % Criterion 3 has been chosen: Maximum deviation from IPS condition
        criterion_name = 'Maximum deviation of IPS condition (difference of middle valued eigenvalue to one)';
        if handles.asc_status_ILS(3) == 0
            handles.asc_number_ILS = handles.asc_number_ILS + 1; % increase number of asc
            default_value = 1.e-2;
            handles.asc_list_ILS(handles.asc_number_ILS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_ILS(3) = handles.asc_number_ILS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_ILS, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
        %
    case 4 % Criterion 4: Maximum misorientation of CPPs {110}_alpha and {111}_gamma
        criterion_name = 'Maximum rotation angle of lath inclusion';
        if handles.asc_status_ILS(4) == 0
            handles.asc_number_ILS = handles.asc_number_ILS + 1; % increase number of asc
            default_value = 10.0;
            handles.asc_list_ILS(handles.asc_number_ILS) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_ILS(4) = handles.asc_number_ILS; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_ILS, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
end

% --- Executes on selection change in popup_asc_block_level.
function popup_asc_block_level_Callback(hObject, eventdata, handles)

selected_crit = get(hObject,'Value');
switch selected_crit
    case 1 % Criterion 1 : Minimum slip plane density
        criterion_name = 'Rotation angle of (average) block deformation.';
        if handles.asc_status_blocks(1) == 0 % if inactive
            handles.asc_number_blocks = handles.asc_number_blocks + 1; % increase number of asc
            default_value = 2.0;
            handles.asc_list_blocks(handles.asc_number_blocks) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_blocks(1) = handles.asc_number_blocks; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_blocks, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
        %
    case 2 % Criterion 1 : Minimum slip plane density
        criterion_name = 'Deviation of IPS condition (middle valued eigenvalue to one).';
        if handles.asc_status_blocks(2) == 0 % if inactive
            handles.asc_number_blocks = handles.asc_number_blocks + 1; % increase number of asc
            default_value = 1.e-2;
            handles.asc_list_blocks(handles.asc_number_blocks) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_blocks(2) = handles.asc_number_blocks; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_blocks, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
        %
    case 3 % Criterion 1 : Minimum slip plane density
        criterion_name = 'Deviation of block habit plane to {111}_aust';
        if handles.asc_status_blocks(3) == 0 % if inactive
            handles.asc_number_blocks = handles.asc_number_blocks + 1; % increase number of asc
            default_value = 5.0;
            handles.asc_list_blocks(handles.asc_number_blocks) = selected_crit; % keep track of which criterion is at which point in the asc list
            handles.asc_status_blocks(3) = handles.asc_number_blocks; % set = number in row, in order to show that crit is already active and when it is to be applied
            handles = create_asc_panel_MartCalc(handles, handles.pan_asc_blocks, criterion_name, default_value, selected_crit );
            guidata(hObject, handles); % Update handles structure
        else
            updateLog_MartCalc(hObject, handles, ['Criterion - "',criterion_name,'" is already active!'])
        end
        %
end


%#########################################################################
%% --- Executes on button press in pushbutton_update_lath_crit.
function pushbutton_update_lath_crit_Callback(hObject, eventdata, handles)
%
set(handles.InterfaceObj,'Enable','off');
update_lath_selection_criteria;
set(handles.InterfaceObj,'Enable','on'); 

% --- Executes on button press in pushbutton_update_block_crit.
function pushbutton_update_block_crit_Callback(hObject, eventdata, handles)

set(handles.InterfaceObj,'Enable','off');
update_block_selection_criteria;
set(handles.InterfaceObj,'Enable','on'); 


%#########################################################################
%% --- Executes on selection change in popup_sorting_IPS.
function popup_sorting_IPS_Callback(hObject, eventdata, handles)
%
%
% reduced_solutions are now directly availble since volumetric changes are
% thrown out immediately. -> reduced solutions
if isempty(handles.martensite.IPS_solutions.array)
    updateLog_MartCalc(hObject, handles,'No solutions available for sorting.')
else
    % disable interface during function call
    set(handles.InterfaceObj,'Enable','off');
    unsrt_sols = handles.reduced_solutions_IPS;
    selected_crit = get(hObject,'Value');
    try
        switch selected_crit
            case 1
                handles.reduced_solutions_IPS = unsrt_sols.sort( 'stepwidth' );
            case 2
                handles.reduced_solutions_IPS = unsrt_sols.sort( 'eps_ips' );
            case 3
                handles.reduced_solutions_IPS = unsrt_sols.sort( 'theta_CPPs' );
            case 4
                handles.reduced_solutions_IPS = unsrt_sols.sort( 'theta_h_to_CPP' );
            case 5
                handles.reduced_solutions_IPS = unsrt_sols.sort( 'theta_KS_min' );
            case 6
                handles.reduced_solutions_IPS = unsrt_sols.sort( 'theta_NW_min' );
            case 7
                handles.reduced_solutions_IPS = unsrt_sols.sort('theta_max_ILSdir_to_h');
                %         case 8
                %             handles.reduced_solutions_IPS = unsrt_sols.sort( 'delta_determinant_max' );
        end
        updateLog_MartCalc(hObject, handles,'Sorting finished.')
    catch ME
        updateLog_MartCalc(hObject, handles,ME.message)
    end
    set(handles.InterfaceObj,'Enable','on');
end
guidata(hObject, handles);
% enable interface again



% --- Executes on selection change in popup_sorting_ILS.
function popup_sorting_ILS_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sorting_ILS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sorting_ILS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sorting_ILS



%#########################################################################
%% BLOCK PART %%

% --- Executes on button press in pushbutton_start_block_calc.
function pushbutton_start_block_calc_Callback(hObject, eventdata, handles)

% disable interface during function call
%set(handles.InterfaceObj,'Enable','off');
updateLog_MartCalc(hObject, handles, '------------- Retrieving input from GUI --------------');
% read user input from GUI for determination of solutions

% Minors relation tolerances
det_tol = num2str(handles.edit_minors_det.String);  
cof_tol = num2str(handles.edit_minors_cof.String);
    
switch handles.popup_calc_lath_level.Value
    case 1 % IPS laths
        %
        if isempty(handles.martensite.IPS_solutions.array) 
            updateLog_MartCalc(hObject, handles, 'Lath IPS-solutions must be calculated first')
        else
            calculation_method = 'Checking minros relations and build bi-variant blocks from lath IPS-solutions (+ optimized phase fractions)';
            updateLog_MartCalc(hObject, handles, [calculation_method,' - started, calculating...']);
            %
            handles.martensite.block_solutions_from_IPS =  mixing_of_atomic_level_solutions( handles.reduced_solutions_IPS, handles.martensite.U, cof_tol, det_tol);
            %
            updateLog_MartCalc(hObject, handles, ['Determination of bi-variant blocks fullfilling minors relations from lath IPS solutions completed.', ...
                num2str(length(handles.martensite.block_solutions_from_IPS)),' possibilites within tolerance range.'] );
        end
    case 2 % ILS laths
        %
        if isempty(handles.martensite.ILS_solutions.array) 
            updateLog_MartCalc(hObject, handles, 'Lath ILS-solutions must be calculated first')
        else
            calculation_method = 'Checking minros relations and build bi-variant blocks from lath ILS-solutions (+ optimized phase fractions)';
            updateLog_MartCalc(hObject, handles, [calculation_method,' - started, calculating...']);
            %
            handles.martensite.block_solutions_from_ILS =  mixing_of_atomic_level_solutions( handles.reduced_solutions_ILS, handles.martensite.U, cof_tol, det_tol);
            updateLog_MartCalc(hObject, handles, ['Determination of bi-variant blocks fullfilling minors relations from lath ILS solutions completed.', ...
                num2str(length(handles.martensite.block_solutions_from_ILS)),' possibilites within tolerance range.'] );
        end
        
        %%            Qi 2014 Acta - removed from direct block calculation
        %             get_input_MartCalc;
        %             calculation_method = 'direct block approach, mirrorsym. & equal double-shears';
        %             updateLog_MartCalc(hObject, handles, [calculation_method,' - started, calculating...']);
        %             %updateLog_MartCalc(hObject, handles, 'please wait...');
        %             %
        %             % highly symmetric mirror planes from bcc
        %             % {001} family
        %             sort_out_negatives = true;
        %             ms = all_from_family_perms( [0 0 1], sort_out_negatives );
        %             % {011} family
        %             ms = cat(1, ms, all_from_family_perms( [0 1 1], sort_out_negatives ) );
        %             handles.martensite.mirror_planes = ms;
        %             %
        %             handles.martensite.IPS_solutions = block_symmetric_doubleshear(handles.martensite, handles.austenite);
        %             updateLog_MartCalc(hObject, handles, ['Determination of (direct) composite block solutions completed: ' num2str(size(handles.martensite.IPS_solutions.array,2)),' solutions found.'] );
        %             %
        %             update_lath_selection_criteria;
end

guidata(hObject, handles);
% enable interface again
%set(handles.InterfaceObj,'Enable','on');

   

%#########################################################################
%% --- Executes on button press in pushbutton_write_solutions.
function pushbutton_write_solutions_Callback(hObject, eventdata, handles)

if isfield(handles,'reduced_solutions_IPS')
    % write lath IPS solutions
    if isempty(handles.reduced_solutions_IPS.array) % ~handles.reduced_solutions.solutions_available
        updateLog_MartCalc(hObject, handles, 'No lath solutions available.');
    else
        filename = handles.filename_results_edittext.String; %{1};
        write_input_parameters(filename,'w', handles.martensite, handles.austenite);
        %
%        if isfield(handles,'reduced_solutions')
            write_calc_specs_laths(filename, 'a',     handles.martensite, handles.reduced_solutions_IPS);
            write_lath_solutions(filename, 'a', handles.martensite, handles.reduced_solutions_IPS);
%        else
%            write_calc_specs_laths(filename, 'a', handles.martensite);
%            write_lath_solutions(filename, 'a', handles.martensite);
%        end
    end % end if solutions_available = true
    
    % write optimized block solutions from lath solutions
    if isfield(handles,'block_solutions')
        write_opt_block_solutions(filename, 'a',handles.block_solutions)
    end
    updateLog_MartCalc(hObject, handles, 'Writing results finished.');
    
else
    updateLog_MartCalc(hObject, handles, 'Please execute lath calculations and reduction first!');
end




%#########################################################################
%% check input of EDIT TEXT fields directly after typing it in

function lc_edtxt_aust_val_Callback(hObject, eventdata, handles)
% if input is a string but not a number change it back to the default value
%if isnan(str2double(get(hObject,'String') ) )
%    set(hObject,'String','3.6017264');   
%end
if ~str2double( get( handles.lc_edtxt_aust_val,'String') ) > 0
    set(hObject,'String','3.6017264');
    updateLog_MartCalc(hObject, handles,'No reasonable input for austenite lattice parameter - please correct!');  
end



function lc_edtxt_mart_val_Callback(hObject, eventdata, handles)
% if input is a string but not a number change it back to the default value
%if isnan(str2double(get(hObject,'String') ) )
%  set(hObject,'String','2.8807346');
%end
if ~str2double( get( handles.lc_edtxt_mart_val,'String') ) > 0
    set(hObject,'String','2.8807346');
    updateLog_MartCalc(hObject, handles,'No reasonable input for austenite lattice parameter - please correct!');  
end


function edit_minors_det_Callback(hObject, eventdata, handles)
% 
det_tol = str2double(get(hObject,'String') );
if isnan(det_tol)
    set(hObject,'String','1.e-3');
else
    if( det_tol > 1.e-3 )
        updateLog_MartCalc( hObject, handles,'A higher tolerance than 1.e-3 for minors relations is not allowed.');
        set(hObject,'String','1.e-3');
    end
end

function edit_minors_cof_Callback(hObject, eventdata, handles)
%  
cof_tol= str2double(get(hObject,'String') );
if isnan(cof_tol)
    set(hObject,'String','1.e-3');
else
    if( cof_tol > 1.e-3 )
        updateLog_MartCalc( hObject, handles,'A higher tolerance than 1.e-3 for minors relations is not allowed.');
        set(hObject,'String','1.e-3');
    end
end

% function max_misori_HPs_laths_for_blocks_edit_txt_Callback(hObject, eventdata, handles)
% % if input is a string but not a number change it back to the default value
% in = str2double(get(hObject,'String') );
% if isnan( in ) || in > 90
%   set(hObject,'String','');
% end
% 
% function misori_HPintersec_cpdir_edit_txt_Callback(hObject, eventdata, handles)
% % if input is a string but not a number change it back to the default value
% in = str2double(get(hObject,'String') );
% if isnan( in ) || in > 90
%   set(hObject,'String','');
% end


function filename_results_edittext_Callback(hObject, eventdata, handles)
% bla
filename = get(hObject,'String');
%if ~isempty(regexp(fname, ['^(?!^(PRN|AUX|CLOCK\$|NUL|CON|COM\d|LPT\d|\..*)', ...
%        '(\..+)?$)[^\x00-\x1f\\?*:\"><|/]+$'], 'once') )
if ~isempty(regexp(filename, '[/\*:?"<>|]', 'once'))
    set(hObject,'String','');
end

set(handles.InterfaceObj,'Enable','on');

%handles.reduced_solutions_ILS.array.u
%handles.reduced_solutions_ILS.array(5).slip.eps_s
%handles.reduced_solutions_ILS.array.slip.slip_normal_plane_vec
%handles.reduced_solutions_ILS.array.slip.shear_direction

% save('ILS_solutions_line_search_delta_eps_ini_0_3__min_delta_res_1_e-3', 'handles.reduced_solutions_ILS' ) ;
% handles.reduced_solutions_ILS.array.shear_increments









%#########################################################################
%% all functions after this are not used...

% --- Executes on selection change in mixing_criteria_for_blocks.
function mixing_criteria_for_blocks_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function lsc_popup_IPS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lsc_popup_IPS (see GCBO)
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


% --- Executes during object creation, after setting all properties.
function popup_calc_lath_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_calc_lath_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popup_asc_block_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_asc_block_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function mixing_criteria_for_blocks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mixing_criteria_for_blocks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popup_sorting_IPS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sorting_IPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function filename_results_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename_results_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit156_Callback(hObject, eventdata, handles)
% hObject    handle to filename_results_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename_results_edittext as text
%        str2double(get(hObject,'String')) returns contents of filename_results_edittext as a double


% --- Executes during object creation, after setting all properties.
function edit156_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename_results_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit157_Callback(hObject, eventdata, handles)
% hObject    handle to edit157 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit157 as text
%        str2double(get(hObject,'String')) returns contents of edit157 as a double





% --- Executes during object creation, after setting all properties.
function max_misori_HPs_laths_for_blocks_edit_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_misori_HPs_laths_for_blocks_edit_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit159_Callback(hObject, eventdata, handles)
% hObject    handle to edit159 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit159 as text
%        str2double(get(hObject,'String')) returns contents of edit159 as a double


% --- Executes during object creation, after setting all properties.
function edit159_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit159 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function misori_HPintersec_cpdir_edit_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to misori_HPintersec_cpdir_edit_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function lsc_popup_ILS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lsc_popup_ILS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit161_Callback(hObject, eventdata, handles)
% hObject    handle to edit161 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit161 as text
%        str2double(get(hObject,'String')) returns contents of edit161 as a double


% --- Executes during object creation, after setting all properties.
function edit161_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit161 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit162_Callback(hObject, eventdata, handles)
% hObject    handle to edit162 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit162 as text
%        str2double(get(hObject,'String')) returns contents of edit162 as a double


% --- Executes during object creation, after setting all properties.
function edit162_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit162 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit163_Callback(hObject, eventdata, handles)
% hObject    handle to edit163 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit163 as text
%        str2double(get(hObject,'String')) returns contents of edit163 as a double


% --- Executes during object creation, after setting all properties.
function edit163_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit163 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit164_Callback(hObject, eventdata, handles)
% hObject    handle to edit164 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit164 as text
%        str2double(get(hObject,'String')) returns contents of edit164 as a double


% --- Executes during object creation, after setting all properties.
function edit164_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit164 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function popup_sorting_ILS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sorting_ILS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit169_Callback(hObject, eventdata, handles)
% hObject    handle to edit169 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit169 as text
%        str2double(get(hObject,'String')) returns contents of edit169 as a double


% --- Executes during object creation, after setting all properties.
function edit169_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit169 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


    
% --- Executes during object creation, after setting all properties.
function edit_minors_det_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minors_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit171_Callback(hObject, eventdata, handles)
% hObject    handle to edit171 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit171 as text
%        str2double(get(hObject,'String')) returns contents of edit171 as a double


% --- Executes during object creation, after setting all properties.
function edit171_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit171 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_minors_cof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minors_cof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
