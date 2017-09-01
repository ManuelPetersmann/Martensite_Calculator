function[reduced_planes, reduced_directions] = reduce_ambiguous_slipsystems(all_planes, all_directions, count_directions_extra)
% takes Nx3 arrays "all_planes" and "all_directions" a list of slip planes and directions 
% defining slip systems { directions(i,:) , planes(i,:) }
% and sorts out ambiguously occuring systems. 
% the third argument is a boolen to specify wheter each deformation
% direction of a slip system should be counted extra, which may be useful
% in calculations, however not to determine e.g. the number of independent
% slip systems

if nargin < 3
    count_directions_extra = false;
end

current_systems = 0;
reduced_planes = ones(1,3);
reduced_directions = ones(1,3);

for x = 1:size(all_planes,1)
    contained = false;
    m_e = all_planes(x,:);
    l_e = all_directions(x,:);
    
    % if its a valid system - check if this system has already been found
    for y = 1:size(reduced_planes,1)
        m_a = reduced_planes(y,:);
        l_a = reduced_directions(y,:);
        % sort out ambiguous combinations (4 possibilities for one system
        % because of opposite directions)
        % cross product to determine same but differently oriented ones
        % yield same results as checking the differences
        %if ( abs( cross(m_e,m_a) ) < 1e-4) & ( abs( cross(l_e,l_a) ) < 1e-4 )
        
        if ( count_directions_extra && ( ( sum(abs(m_e - m_a)) < 1e-5 ) || ( sum(abs(m_e + m_a)) < 1e-5 ) ) )
            if (sum(abs(l_e - l_a)) < 1e-5 ) 
            contained = true;
            break  
            end
        elseif ( ( sum(abs(m_e - m_a)) < 1e-5 ) || ( sum(abs(m_e + m_a)) < 1e-5 ) )
            if ( (sum(abs(l_e - l_a)) < 1e-5 ) || ( sum(abs(l_e + l_a)) < 1e-5 ) )
                contained = true;
                break
            end
        end

          
    end % end for
    
    % if it has not been found yet add it to independent oness
    if ~ contained
        reduced_planes(current_systems+1,:) = m_e; 
        reduced_directions(current_systems+1,:) = l_e;
        current_systems = current_systems + 1;
    end
    
end % for end
end % function end

