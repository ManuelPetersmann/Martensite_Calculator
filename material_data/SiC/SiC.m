clear all
clc

%% Set up coordinate systems and lattice parameters
martensite = Martensite(); % creates martensite object
austenite = Base();

% conventional bases
%austenite.my_base =  [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
%martensite.my_base = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
austenite.Bravais_type  = 'cubic'; % fcc
martensite.Bravais_type = 'hexagonal'; % hcp
% austenite.Centering = 
% martensite.Centering = 'I';  % space group is: P63mmc

x_sample = [1 0 0]';
y_sample = [0 1 0]';
z_sample = [0 0 1]';

%alpha = [ 210 ];
%beta =  [ 102 ];
%gamma = [ 77 ];

tolerance_lambda_2 = 0.005

%%

% orientation of dominate habit plane in paper plane specimen Phase_5_Zentrum_etc
HPangle = deg2rad(80); % 80 deg measured from image 
min_delta = 1.; % in radians - variable for minimum


%% Lattice Correspondence Matrix and structural stretch tensor from structtrans (www.structrans.org):
% e_mart = C_am * e_aust; infinetively many choices. Take one for slip system transformation....

% this corresponds only to U1 !!!
martensite.C_am = [0.5, 0.5, 0;    % [ 1 0 0 ]
                   0.,  0.5, 0.5;  % [ 0 1 0 ]
                   1.5, -1., 1.5]; % [ 0 0 1 ]
%%correspondance_matrix(austenite, martensite, m1, m2);

U1 = [	 -0.0383,	-0.0267,	-0.0374;
   	  	-0.0267,	 1.0840,	-0.0267;
	  	-0.0374,	-0.0267,	 0.9617 ];

martensite.U = U1;

Us = martensite.symmetry_related( martensite.U )
Us = reduce_matrices( Us )

display(['Volume change in percent is:', num2str( ( det(U1)-1)*100 ) ] );
%vars = martensite.variants()

dirstrain = directional_strain( [1., 1., 1.,], U1 );
display(['strain in [1-11]_fcc habit plane direction is ', num2str( dirstrain )] );

dirstrain = directional_strain( [1., 1., 0.,], U1 );
display(['strain in [110]_fcc direction normal to habit plane is ', num2str( dirstrain )] );

[E,V] = eigs(U1)

%% other Stretch tensor solution (see data!)
% U2 = [];
% display(['Volume change of second most promising U in percent is:', num2str( ( det(U2)-1)*100 ) ] );
%cp = B3*C_am; 
%%

%[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( U1'*U1 )
%tolerance = 1.e-4;
% [ is_possible_solution , lambda2_smaller1(is) ] = check_IPS_solution(lambda_1, lambda_2, lambda_3, tolerance);

%farbNr = 0;
%farben = ['c','r','g','b','y','m']; % ['r','g','b','c','m', 'y', 'k']; 
C = {'c','m','k','b','r','g','y',[.5 .6 .7],[.8 .2 .6], [.1 .9 .5], [.3 .3 .9], [.9 .1 .7]};
% zeichen = ['o', 's', 'x', '.', 'd', '^', 'p', 'h'];
lines = ['-', '--']; % , ':', ':','-.'];

figure;
isol = 0;
k = 0; % counter for plotting colors
% 1 - lambda_2 = 0.0048  - im Vergleich zu TiAl publication 0.04 - also fast factor 10

array_index = 1;
for i = 1:size(Us,3) 
    
    k = k+1; % for line color
    
    for j = 1: length(alpha)

    rot_mat(:,:,i) = eulerAngleMatrix( alpha(j), beta(j), gamma(j), 'bungee');
    
    % rotates sample directions into crystal directions by:
    % e_crystal_i = rot_mat * e_sample_i
    % the direction that we want to be transformed is [001] sample.
    
    % x = A\b  % Denotes the solution to the matrix equation Ax = b, obtained using mldivide.
    % system of 9 linear equations    rot_mat' * e_crystal = e_sample
    
    A = rot_mat(:,:,i)';
        
    
    %% actually there should be 4 hps!
    %[y1, y3, d1, d2, h1, h2, Q1, Q2] = rank_one( U1 , eye(3), tolerance_lambda_2 );
    [y1, y3, d1(i,:), d2(i,:), h1(i,:), h2(i,:), Q1(:,:,i), Q2(:,:,i)] = rank_one( Us(:,:,i) , eye(3), tolerance_lambda_2 );
%    d1(i,:)
%    d2(i,:)
    if i == 1
        h1(i,:)
        h2(i,:)
    end
    eps_0 = sqrt(y3) - sqrt(y1);
    R = Q1(:,:,i);
    ang = acos( ( trace(Q1(:,:,i) ) -1.) / 2.); 
    axis = 1./(2*sin(ang));
    axis = axis * [R(3,2)-R(2,3)  ,  R(1,3)-R(3,1)  ,  R(2,1)-R(1,2)]; 
    ang = acos( ( trace(Q2(:,:,i) ) -1.) / 2.);
    R = Q2(:,:,i);
    axis = 1./(2*sin(ang));
    axis = axis * [R(3,2)-R(2,3)  ,  R(1,3)-R(3,1)  ,  R(2,1)-R(1,2)]; 
    
    isol = isol + 2; % increase counter for number of solutions found
    % h1(i,:)
    % h2(i,:)
    % Create Slip_solution objects and append them to object array
    %solutions.array( isol-1 ) =  Slip_solution_multishear(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, shear_magsfound, dsfound, nsfound );
    %solutions.array( isol )   =  Slip_solution_multishear(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, shear_magsfound, dsfound, nsfound );
    % solutions.array( isol-1 ) =  Slip_solution(F, I, y1, y3, d1, h1, Q1, Q1*martensite.U, eps_s, d, n );
    % solutions.array( isol )   =  Slip_solution(F, I, y1, y3, d2, h2, Q2, Q2*martensite.U, eps_s, d ,n );
    % solutions.array( isol-1 ).id = isol-1;
    % solutions.array( isol ).id = isol;

%% calculation in crystal coordinate system (more difficult - checked that it is the same! )    
%     x_p(:,i) = A \ x_sample;
%     y_p(:,i) = A \ y_sample;
%     z_p(:,i) = A \ z_sample;
%     % frame(:,:,i) = cat( 2,x_p(:,i),y_p(:,i),z_p(:,i) );

%     schnitt3D(i,1:3) = cross(h1(i,:), z_p(:,i)');
%     schnitt3D(i,4:6) = cross(h2(i,:), z_p(:,i)'); %z_p - paper plane
%     schnitt2D(i,1) = dot(x_p(:,i), schnitt3D(i,1:3) );
%     schnitt2D(i,2) = dot(y_p(:,i), schnitt3D(i,1:3) );
%     schnitt2D(i,3) = dot(x_p(:,i), schnitt3D(i,4:6) );
%     schnitt2D(i,4) = dot(y_p(:,i), schnitt3D(i,4:6) );
%     schnitt2D(i,1:2) = schnitt2D(i,1:2) / norm( schnitt2D(i,1:2) );
%     schnitt2D(i,3:4) = schnitt2D(i,3:4) / norm( schnitt2D(i,3:4) );
%     
%     theta1 = atan2(schnitt2D(i,1),schnitt2D(i,2));
%     theta2 = atan2(schnitt2D(i,3),schnitt2D(i,4));
%     
%     myPolar([theta1,theta1],[-1,1] ,1, C{i} );
%     hold all
%     myPolar([theta2,theta2],[-1,1] ,1, C{i} );
%     figure;

%          if mod(i,orientations_per_grain) == 0
%              k = k + 1
%          end
    
    h1_sample = A * h1(i,:)';
    h2_sample = A * h2(i,:)';
    d1_sample = A * d1(i,:)';
    d2_sample = A * d2(i,:)';
    F1_sample = eye(3) + eps_0 * d1_sample' * h1_sample;  
    F2_sample = eye(3) + eps_0 * d2_sample' * h2_sample; 
    e1_sample = 0.5 * (F1_sample'*F1_sample) - eye(3);
    e2_sample = 0.5 * (F2_sample'*F2_sample) - eye(3);
    
    %% calculation in sample coordinate system
    % cut with polished specimen surface
    th1 = cross(h1_sample, [0. 0. 1.] );
    th2 = cross(h2_sample, [0. 0. 1.] );
    % PET 24.05.2018 -- atan2( Y,X )  --- corrected Y,X as in help atan2  (exchanged order )
    % yields values that can be compared to the variable "angle" in rad
    th11 = atan2( th1(2), th1(1) );
    th22 = atan2( th2(2), th2(1) );
    delta11 = abs(HPangle - th11);
    delta22 = abs(HPangle - th22);
    
    deltas(array_index) = delta11;
    eners(array_index) = e1_sample(3,3);
    array_index = array_index + 1;
    deltas(array_index) = delta22;
    eners(array_index) = e2_sample(3,3);
    array_index = array_index + 1;
    
    if delta11 > delta22
        delta = delta22;
    else
        delta = delta11;
    end
    if min_delta > delta
        min_delta = delta;
    end
    myPolar([th11,th11],[-1,1] ,1, C{k} );
    hold all
    myPolar([th22,th22],[-1,1] ,1, C{k} );
    %figure;
    
    if i == 1
        view(-90,90); % rotate plot such that it is aligned like EBSD image
        %hold all
    end
%% END Script
    
%    farbNr = farbNr +1;
%    farbe = farben( farbNr );
%     if farbNr == 4
%         farbNr = 0
%     end
    
    %     if mod(i,2*shears) == 1
    %         farbNr = farbNr +1;
    %         farbe = farben( farbNr );
    %

    %         figure('position', [0, 0, 500, 500])
    %         axis ([-1, 1, -1, 1])
    %         axis square
    %         angles = linspace(0,2*pi);
    %         plot( cos(angles), sin(angles),'b');
    %         view(90,90) % rotate the plot clockwise by 90 degree to coincide with the EBSD coordinate system
    %         hold all
    %    end
    
%     if mod(i, 2) == 0
%         line = lines(2);
%     else
%         line = lines(1);
%     end    
    
%    einstellungen = [rand(1,3), line]; %line, farbe, zeichen]
    
    %     if i == 1
    %         figure('position', [0, 0, 500, 500])
    %         axis ([-1, 1, -1, 1])
    %         axis square
    %         angles = linspace(0,2*pi);
    %         plot( cos(angles), sin(angles),'b');
    %         view(90,90) % rotate the plot clockwise by 90 degree to coincide with the EBSD coordinate system
    %         hold all
    %     end
    
    %plot([xy(i,1), -xy(i,1)],[xy(i,2), -xy(i,2) ], einstellungen, 'MarkerSize',5, 'LineWidth',2) %  'MarkerSize',4,'MarkerFaceColor','b' 'LineWidth',1, ,  'LineWidth',1


% rotate plot
%view(90,90)
    end

end

deltas_eners = cat(1,deltas, eners);
deltas_eners'

