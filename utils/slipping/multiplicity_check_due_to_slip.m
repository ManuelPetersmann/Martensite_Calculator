function sols = multiplicity_check_due_to_slip(lath_solutions)
% call: multiplicity_check_due_to_slip(lath_solutions)

% loops over .F1 property and checks if same matrices appear.
% if so - add shear - slip information (slip systems and shear magnitudes)
% to first solution where F1 has been found.

multiplicity = 1;
sols = lath_solutions.array;
start_length = size(lath_solutions.array,2);

% reduced_by = 0;

% loop over slip system combinations
% for is1 = 1: (size(lath_solutions.array,2)-1)
%     for is2 = (is1+1): size(lath_solutions.array,2) % for loop index
%     cannot be set withing the for loop - set only once...
is1 = 0;
while is1 < size(sols,2)-1
    is1 = is1 + 1;
    F1 = sols(is1).F1;
    is2 = is1 + 1;
    while is2 < size(sols,2)
        is2 = is2 + 1;

        F2 = sols(is2).F1;
        % Considering that the two F could be equal since the slip
        % deformations are not linearly independent! c.f. non-uniqueness of
        % plastic slip
        if sum(sum(abs(F1 - F2 ))) < 1.e-3
            % safe slip information from sol2 to sol1
            multiplicity = multiplicity + 1;
            sols(is1).id(multiplicity) = sols(is2).id;
            sols(is1).eps_s(:,:,multiplicity) = sols(is2).eps_s;
            sols(is1).shear_direction(:,:,multiplicity) = sols(is2).shear_direction;
            sols(is1).slip_normal_plane_vec(:,:,multiplicity) = sols(is2).slip_normal_plane_vec;

            % remove multiply occuring solution from array
            sols(is2) = []; 
%             reduced_by = reduced_by + 1;
        end
        
    end % end of loop 1
end % end of loop 2
end_length = size(sols,2);

disp( ['Solutions reduced from ',num2str(start_length),' to ',num2str(end_length),' due to equal deformations, only different active slip -> non-uniqueness, linear dependence of slip deformation'] )

end
