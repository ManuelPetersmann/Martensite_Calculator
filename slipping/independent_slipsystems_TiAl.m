function [ m_red, l_red ] = independent_slipsystems( R )
% determines the number of slip systems ( without ambiguities ) in a 
% system derived from variant one


%% slip systems b.c.c or h.c.p
% since the shear is a substantial part of the transformation only some 
% CORRESPONDING shear systems which are favorable in the cubic as well as 
% in the hexagonal lattice are considered. These are:

l = [1 1 1]; % #1
m = [-1 1 0];
l(2,:) = [1 1 1]; % #2
m(2,:) = [0 1 -1];
l(3,:) = [1 1 1]; % #3
m(3,:) = [1 0 -1];
l(4,:) = [1 1 1]; % #4
m(4,:) = [1 1 -2];
l(5,:) = [1 1 1]; % #5
m(5,:) = [1 -2 1];
l(6,:) = [1 1 1]; % #6
m(6,:) = [-2 1 1];
l(7,:) = [1 -1 1]; % #7
m(7,:) = [0 1 1];
l(8,:) = [1 -1 1]; % #8
m(8,:) = [1 0 -1];
l(9,:) = [1 -1 1]; % #9
m(9,:) = [1 1 0];
l(10,:) = [1 -1 1]; % #10
m(10,:) = [1 2 1];
l(11,:) = [1 -1 1]; % #11
m(11,:) = [-1 1 2];
l(12,:) = [0 1 -1]; % #12
m(12,:) = [0 1 1];
l(13,:) = [0 1 -1]; % #13
m(13,:) = [1 1 1];
l(14,:) = [0 1 -1]; % #14 
m(14,:) = [2 1 1];
l(15,:) = [0 1 -1]; % #15
m(15,:) = [3 1 1];
l(16,:) = [0 1 0]; % #16 explicit solution cannot be found while solving
% the quadratic equation for g
m(16,:) = [1 0 1];
l(17,:) = [0 1 0]; % #17
m(17,:) = [0 0 1];
l(18,:) = [0 1 0]; % #18
m(18,:) = [1 0 0];
l(19,:) = [0 1 0]; % #19
m(19,:) = [1 0 2];
l(20,:) = [1 0 0]; % #20
m(20,:) = [0 1 1];
l(21,:) = [1 0 0]; % #21
m(21,:) = [0 1 0];
l(22,:) = [1 0 0]; % #22
m(22,:) = [0 0 1];
l(23,:) = [1 0 0]; % #23
m(23,:) = [0 1 2];
l(24,:) = [1 1 3]; % #24
m(24,:) = [1 -1 0];
l(25,:) = [1 1 3]; % #25
m(25,:) = [1 2 -1];
l(26,:) = [-1 1 3]; % #26
m(26,:) = [1 1 0];
l(27,:) = [-1 1 3]; % #27
m(27,:) = [2 -1 1];
l(28,:) = [-1 1 3]; % #28
m(28,:) = [1 -2 1];
l(29,:) = [3 -1 1]; % #29
m(29,:) = [0 1 1];
l(30,:) = [3 -1 1]; % #30
m(30,:) = [1 1 -2];
l(31,:) = [3 -1 1]; % #31
m(31,:) = [1 2 -1];
% Ergänztes schersystem nach Otte (dieses system hat nichts an den
% Ergebnissen geändert.
l(32,:) = [-1 1 0]; % #32
m(32,:) = [1 1 -2];


% m = [-1 0 1; 1 1 0]; % = {-1 0 1}
% l = [1 1  1; 1 1 1]; % = < 1 1 1 >

mn = [0 0 0];
ln = [0 0 0];
k = 0;
for i = 1:size(m,1)
    for j = 1:size(R,3)
        k=k+1;
        % Drehung des Vektors bei festem Koordinatensystem
        mn(k,:) = R(:,:,j) * m(i,:)'; % Drehung des Koordinatensystems für Vektor wäre R^-1*v
        ln(k,:) = R(:,:,j) * l(i,:)'; 
    end
end

%% or initialisation via permutations of expected indices of miller index families
% m1 = [-1 0 1]; % = {-1 0 1}
% m2 = [1 1 0];
% l1 = [1 1 1]; % = < 1 1 1 >
% %l2 = [1 -1 1];
% 
% 
% mm = cat(1, perms( m1 ), perms( m2 ) );
% ll = perms( l1 ); %, perms( l2 ) );
% 
% nr = 1;
% for i = 1:size(mm,1)
%     for j = 1:size(ll,1)
%         mn(nr,:) = mm(i,:); 
%         ln(nr,:) = ll(j,:);
%         nr = nr +1;
%     end
% end

%% In any case a reduction to unequivalent systems has to be made
[m_red, l_red] = reduce_ambiguous_slipsystems( mn, ln );
length( m_red )
length( l_red )
% m_red = mn
% l_red = ln

end


