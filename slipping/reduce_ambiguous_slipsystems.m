function[reduced_planes, reduced_directions] = reduce_ambiguous_slipsystems(all_planes, all_directions)
% takes a list of slip planes and directions each as a Nx3 Matrix
% defining slip systems { directions(i,:) , planes(i,:) }
% and sorts out ambiguously occuring systems

current_systems = 0;
reduced_planes = ones(1,3);
reduced_directions = ones(1,3);

for x = 1:size(all_planes,1)
    contained = false;
    m_e = all_planes(x,:);
    l_e = all_directions(x,:);
    % sort out permutations where m and l are not mutually perpendicular
    % hence this combination is not a valid slip system
    if abs( dot(m_e, l_e) ) > 1e-5
        continue
    end
    
    % if its a valid system - check if this system has already been found
    for y = 1:size(reduced_planes,1)
        m_a = reduced_planes(y,:);
        l_a = reduced_directions(y,:);
        % sort out ambiguous combinations (4 possibilities for one system
        % because of opposite directions)
        % cross product to determine same but differently oriented ones
        % yield same results as checking the differences
        %if ( abs( cross(m_e,m_a) ) < 1e-4) & ( abs( cross(l_e,l_a) ) < 1e-4 )
        if ( ( sum(abs(m_e - m_a)) < 1e-5 ) || ( sum(abs(m_e + m_a)) < 1e-5 ) )
            if ( (sum(abs(l_e - l_a)) < 1e-5 ) || ( sum(abs(l_e + l_a)) < 1e-5 ) )
            contained = true;
            break
            end
        end
    end
    
    % if it has not been found yet add it to independent oness
    if ~ contained
        reduced_planes(current_systems+1,:) = m_e; 
        reduced_directions(current_systems+1,:) = l_e;
        current_systems = current_systems + 1;
    end
    
end % for end
end % function end

