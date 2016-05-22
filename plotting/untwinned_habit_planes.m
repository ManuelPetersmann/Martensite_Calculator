function [ schnitt2D, ns, as, z_p ] = untwinned_habit_planes( variants, eulerAngles, aJm )
% computes possible habit planes given the martensitic variants Ui
% and the orientation of the paper plane. Then plots the cut of the habit
% planes and the paper plane.
% Notion and equations are equal to N. K. Simha: "Twin and Habit Plane Microstructures due
% to the tetragonal to Monoclinic transformation of zirconia

phi1 = (pi/180) * eulerAngles(1);
theta = (pi/180) * eulerAngles(2);
phi2 = (pi/180) * eulerAngles(3);
% get vector normal to paper-plane p
paperplane_coord_syst = eulerAngleRotation( phi1, theta, phi2, [1 0 0; 0 1 0; 0 0 1]);
z_p = paperplane_coord_syst(:,3);
y_p = paperplane_coord_syst(:,2);
x_p = paperplane_coord_syst(:,1);

ns = zeros(1,3);
as = zeros(1,3);

nr = 1;

for i1 = 1:size(variants,3)
        Ui = variants(:,:,i1);
        Ct = Ui^(-2); % wichtig Ui = I setzen! Wegen herleitung von C_ij
        %
        %Next the Eigenvalues and korresponding Eigenvectors of Ct will be calcuated
        [y1, y3, e1, e3] = ordered_eigenvalues(Ct);
        %
        % X = Tschi, r = ro, n = n mit Dach
        %If Ct ~= (not equal) to I, shear a and normal n are given by:
        %
        for i2 = 1:2
            switch i2
                case 1
                    tschi=1; % ACHTUNG GLEICH MINUS zuerst wenn DIREKTER VERGLEICH MIT c++ PROGRAMM DA EIN VEKTOR ENTGEGENGESETZTES VORZEICHEN!!!!
                otherwise 
                    tschi= -1;
            end
            %
            [a,n] = an_bh_calc( y1, y3, e1, e3, tschi, Ui);
            
            % if the shear is formally applied in the parent phase:
            ns(nr,:) = n;
            as(nr,:) = a;
            
            % if the shear is applied in the product phase transform the
            % habit plane normal to the austenite
            
            % mJa = aJm \ eye(3);
            % n = n * mJa ;
            n'
            ns(nr,:) = n';
            % as(nr,:) = ( mJa * a' )';
            
           
            nr = nr +1;           
   
        end
end
        
schnitt3D = zeros(1,3);
schnitt2D = zeros(1,2);
for i = 1:size(ns,1)
    schnitt3D(i,:) = cross(ns(i,:),z_p);
    schnitt2D(i,1) = dot(x_p, schnitt3D(i,:) );
    schnitt2D(i,2) = dot(y_p, schnitt3D(i,:) );
    schnitt2D(i,:) = schnitt2D(i,:) / norm( schnitt2D(i,:) );
end
        
end











