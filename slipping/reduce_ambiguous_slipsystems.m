function[ausgang1, ausgang2] = reduce_ambiguous_slipsystems(eingang1, eingang2)
% takes a list of slip planes and directions each as a Nx3 Matrix
% defining slip systems { directions(i,:) , planes(i,:) }
% and sorts out ambiguously occuring systems

weiter = 1;
ausgang1 = ones(1,3);
ausgang2 = ones(1,3);

for x = 1:size(eingang1,1)
    contained = false;
    m_e = eingang1(x,:);
    l_e = eingang2(x,:);
    % sort out permutations where m and l are not mutually perpendicular
    % hence this combination is not a valid slip system
    if abs( dot(m_e, l_e) ) > 1e-4
        continue
    end
    
    % if its a valid system - check if this system has already been found
    for y = 1:size(ausgang1,1)
        m_a = ausgang1(y,:);
        l_a = ausgang2(y,:);
        % sort out ambiguous combinations (4 possibilities for one system
        % because of opposite directions)
        % cross product to determine same but differently oriented ones
        %if ( abs( cross(m_e,m_a) ) < 1e-4) & ( abs( cross(l_e,l_a) ) < 1e-4 )
        if ( ( sum(abs(m_e - m_a)) < 1e-2 ) | ( sum(abs(m_e + m_a)) < 1e-2 ) )
            if ( (sum(abs(l_e - l_a)) < 1e-2 ) | ( sum(abs(l_e + l_a)) < 1e-2 ) )
            contained = true;
            break
            end
        end
    end
    
    % if it has not been found yet add it to independent oness
    if ~contained
        ausgang1(weiter,:) = m_e; 
        ausgang2(weiter,:) = l_e;
        weiter = weiter + 1;
    end
    
end % for end
end % function end

