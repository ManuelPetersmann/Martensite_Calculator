function[ausgang1, ausgang2] = reduce_ambiguous_slipsystems(eingang1, eingang2)
weiter = 1;
ausgang1 = 678*ones(1,3);
ausgang2 = 678*ones(1,3);
for x = 1:size(eingang1,1)
    u=0;
    m_e = eingang1(x,:);
    l_e = eingang2(x,:);
    % sort out permutations where m and l are not mutually perpendicular
    if abs( dot(m_e, l_e) ) > 1e-4
        continue
    end
    for y = 1:size(ausgang1,1)
        m_a = ausgang1(y,:);
        l_a = ausgang2(y,:);
        % sort out ambiguous combinations (4 possibilities for one system
        % because of opposite directions)
        if ( abs( cross(m_e,m_a) ) < 1e-3) & ( abs( cross(l_e,l_a) ) < 1e-3)
            break
        end
        u = u +1;
    end
    if u == size(ausgang1,1)
        ausgang1(weiter,:) = m_e; % in sort schreiben
        ausgang2(weiter,:) = l_e;
        weiter = weiter + 1;
    end
end % for ende
end % function reduzieren ende

