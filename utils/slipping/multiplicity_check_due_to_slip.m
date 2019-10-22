function lath_solutions = multiplicity_check_due_to_slip(lath_solutions, angle_tol, def_tol) %d1_tol) 
% call: multiplicity_check_due_to_slip(lath_solutions, d1_tol)
%
% Considering that the two F could be equal since the slip
% deformations are not linearly independent! c.f. non-uniqueness of
% plastic slip
% loops over .ST property and checks if same matrices appear.
% the difference between two matrices A and B is measured using the d1 norm:
% sum(sum(abs(A-B)).
% if so - add shear - slip information (slip systems and shear magnitudes)
% to first solution where F1 has been found.
% also sort out solutions with only one active slip system! 
% default d1_tol = 1.e-2
%
% this function is quite slow if the amout of solutions is big (>>1000).
% actually it would be smarter to sort out such solutions with just one
% active slip system (if possible) when assigning solutions to solution
% array...

% if < tolerance then equal (same solution)!
if nargin < 2
%    d1_tol = 1.e-3; % 5.e-3;
    angle_tol = 3.; % 2 degree
    def_tol = 1; %0.5; % percent of relative stretch - for now 1 percent
end

sols = lath_solutions.array;
start_length = size(lath_solutions.array,2);

% reduced_by = 0;

% loop over slip system combinations
% for is1 = 1: (size(lath_solutions.array,2)-1)
%     for is2 = (is1+1): size(lath_solutions.array,2) % for loop index
%     cannot be set withing the for loop - set only once...
is1 = 0;
while is1 < size(sols,2) % <= size(sols,2) -1
    is1 = is1 + 1;
    s1 = sols(is1);
    multiplicity = 1;
    % if one shear magnitude is zero - for now it only works for two slips!
    % check if first shear just has one active system
    zero_shears1 = find( s1.slip.eps_s < 1.e-4 ,1);
    if ~isempty(zero_shears1)
        bigger_zero1 = find(s1.slip.eps_s > 1.e-4,1);
        sl1 = s1.slip.slip_normal_plane_vec( bigger_zero1 , :);
        sd1 = s1.slip.shear_direction(       bigger_zero1 , :);
    end
    %
    is2 = is1 + 1;
    while is2 <= size(sols,2) % < size(sols,2)+1
        s2 = sols(is2);
        
        % if similar / equal within tolerance
        % NOTE that if the F1's are different the ST's could nevertheless
        % be the same, but this case should be distinguished!!!!
%       if sum(sum(abs(s1.F1 - s2.F1 ))) < d1_tol % alternatively ST
%       if sum(sum(abs( inverse(s1.ST)*s2.ST ))) < d1_tol
        nf11 = norm(s1.ST(:,1));
        nf12 = norm(s1.ST(:,2));
        nf13 = norm(s1.ST(:,3));
        nf21 = norm(s2.ST(:,1));
        nf22 = norm(s2.ST(:,2));
        nf23 = norm(s2.ST(:,3));
        if (acosd( dot(s1.ST(:,1) / nf11 , s2.ST(:,1) / nf21 ) ) < angle_tol   && ...
            acosd( dot(s1.ST(:,2) / nf12 , s2.ST(:,2) / nf22 ) ) < angle_tol   && ...
            acosd( dot(s1.ST(:,3) / nf13 , s2.ST(:,3) / nf23 ) ) < angle_tol   && ...
            nf11 - nf21 < (nf11 / 100 )*def_tol   && ...
            nf12 - nf22 < (nf12 / 100 )*def_tol   && ...
            nf13 - nf23 < (nf13 / 100 )*def_tol )   
           
            % check if second shear just has one active system
            zero_shears2 = find( s2.slip.eps_s < 1.e-4, 1);
            
            if ( isempty(zero_shears1)   &&   isempty(zero_shears2) )
                % if matrix difference small and other slip systems, both
                % with slip magnitudes unequal 0!!!
                % safe slip information from sol2 to sol1 then delete sol2
                % sols(is1).id;   sols(is2).id
                multiplicity = multiplicity + 1;
                                
%                 sols(is1).id
%                 s2.id
%                 multiplicity
                sols(is1).id(multiplicity) = s2.id;
                sols(is1).slip.eps_s(:,:,multiplicity) = s2.slip.eps_s;
                sols(is1).slip.shear_direction(:,:,multiplicity) = s2.slip.shear_direction;
                sols(is1).slip.slip_normal_plane_vec(:,:,multiplicity) = s2.slip.slip_normal_plane_vec;
                % also copy ST and
                % LT !!!
                %                
                % remove multiply occuring solution from array
                sols(is2) = [];
                % reduced_by = reduced_by + 1;
            else
                
                % here it is dealt with cases where one or two shears are zero
                if (  ~isempty(zero_shears2)  &&  ~isempty(zero_shears1)  ) % is there a shear with zero magnitude?
                    % in this case entries get only removed
                    bigger_zero2 = find(s2.slip.eps_s > 1.e-4, 1);
                    sl2 = s2.slip.slip_normal_plane_vec(bigger_zero2,:);
                    sd2 = s2.slip.shear_direction(      bigger_zero2,:);
                    % first shear has just one active slip and second has also
                    % just one active one and it is the same --> delete second
                    if cat(2, sl2 == sl1 , sd1 == sd2)
                        sols(is2) = [];
                    end
                end
            end
            
        end % F1-F2 < tol
        
    is2 = is2 + 1;
    
    end % end of loop 1
end % end of loop 2

end_length = size(sols,2);
lath_solutions.array = sols;

disp( ['Solutions reduced from ',num2str(start_length),' to ',num2str(end_length),' due to equal deformations, only different active slip -> non-uniqueness, linear dependence of slip deformation'] )

end
