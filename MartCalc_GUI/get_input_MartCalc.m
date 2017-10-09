% Script for reading input from MartCalc-GUI


%% get lattice parameters
% for austenite lattice
if(num2str(handles.lc_edtxt_aust_val.String) > 0 )
    a_aust = str2num(handles.lc_edtxt_aust_val.String);
else
    updateLog_MartCalc(hObject, handles,'No reasonable input for austenite lattice parameter - please correct!');
    handles.input_status = false;    
    %error('No reasonable input for fcc lattice parameter!');
end
% for martensite lattice
if(num2str(handles.lc_edtxt_mart_val.String) > 0 )
    a_mart = str2num(handles.lc_edtxt_mart_val.String);
else
    updateLog_MartCalc(hObject, handles,'No reasonable input for martensite lattice parameter - please correct!');
    handles.input_status = false; 
    %error('No reasonable input for bcc lattice parameter!');
end
 
%% get input for base vectors from GUI
base_aust = zeros(3,3);
k = 19; % position of 1st entry of 1st base-vec for austenite
for i = 1:3
    for j = 1:3
        base_aust(i,j) = str2num(handles.pan_base_vec.Children(k).String);
        k = k-1;
    end
end
handles.austenite.my_base = base_aust; % checks are done in class!
%
base_mart = zeros(3,3);
k = 9; % position of 1st entry of 1st base-vec for martensite
for i = 1:3
    for j = 1:3
        base_mart(i,j) = str2num(handles.pan_base_vec.Children(k).String);
        k = k-1;
    end
end
handles.martensite.my_base = base_mart; % checks are done in class!

                   
%% get input for correspondence matrix from GUI
C_am = zeros(3,3);
k = 9; % counter for position in array handles.pan_corrmat.Children(k)
for i = 1:3
    for j = 1:3
        C_am(i,j) = str2num(handles.pan_corrmat.Children(k).String);
        k = k-1;
    end
end
handles.martensite.C_am = C_am;

%% Set Bain strain and some other properties
%% EHL: add input in GUI !!! - above all vor CPP and KS (in general set of close packed planes and directions)
handles.austenite.Bravais_type  = 'cubic';
handles.martensite.Bravais_type = 'cubic';
% austenite.Centering = 
handles.martensite.Centering = 'I';

handles.austenite.Lp = a_aust*[1 1 1];  % 3.5975576 % {1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
                                              % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
handles.martensite.Lp = a_mart *[1 1 1];  % 2.8807346  
               
% further initializations
% define Bain-strain
eta1 = (a_mart/a_aust)*sqrt(2);
eta3 = a_mart / a_aust; % this is one form of three pcossible for the bain strain

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf 1.0 getuned
% werden. Die differenz ist also (n1-1).
B3 = [eta1 0    0   
       0  eta1  0
       0  0  eta3];
handles.martensite.U = B3;

updateLog_MartCalc(hObject, handles,['Theoretical volume change of Bain is: ',num2str( (det(handles.martensite.U) -1. )*100),'%'] );

 
%% get input for slip systems
dir_families_aust   = check_input_uitable( handles.uitable_slip_dirs_aust.Data ); % gives 4x3 matrix
plane_families_aust = check_input_uitable( handles.uitable_slip_normals_aust.Data );
dir_families_mart   = check_input_uitable( handles.uitable_slip_dirs_mart.Data );
plane_families_mart = check_input_uitable( handles.uitable_slip_normals_mart.Data );
count_directions_extra = true;
no = zeros(1,3);

%% do verifications on the given input

considered_plasticity = 0;
% if martensite has some input for slip systems
if  ~isequal( dir_families_mart,no ) &&  ~isequal( plane_families_mart, no )
    try
        [handles.martensite.slip_planes, handles.martensite.slip_directions] = independent_slipsystems( plane_families_mart, dir_families_mart, count_directions_extra );
        updateLog_MartCalc(hObject, handles,[num2str(length( handles.martensite.slip_planes )),' shear deformations active in martensite.'] );
        considered_plasticity = 1;
    catch ME
        updateLog_MartCalc(hObject, handles,ME.message);
        handles.input_status = false;
    end
end
% if austenite has some input for slip systems
if  ~isequal( dir_families_aust,no ) &&  ~isequal( plane_families_aust, no )
    try
        [handles.austenite.slip_planes, handles.austenite.slip_directions]   = independent_slipsystems( plane_families_aust,  dir_families_aust,  count_directions_extra );
        updateLog_MartCalc( hObject, handles,[num2str(length( handles.austenite.slip_planes )), ' shear deformations active in austenite.'] );
        considered_plasticity = considered_plasticity + 2;
    catch ME
        updateLog_MartCalc(hObject, handles,ME.message);
        handles.input_status = false;
    end
end

% if neither martensite nor austenite has valid slip systems
if considered_plasticity == 0
    updateLog_MartCalc( hObject, handles,'no valid slip system could be build, please check input' );
    handles.input_status = false;
end

% switch considered_plasticity
% case 1
% case 2
% case 3
% disp( ['Number of possible pairings is = ', num2str( nchoosek(size(ds,1),2) )])
% disp('nr of solutions cannot be greater than 2-times this value.')
% end


%%
% save which phases should be considered for slip
handles.martensite.considered_plasticity = considered_plasticity;

guidata(hObject, handles);


                  
                   
