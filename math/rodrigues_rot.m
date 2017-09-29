% rodrigues_rot - creates 3D-rotation matrix of axis-angle pair (k,theta)
% Direction is determined by the right-hand (screw) rule.
% only possible for for the standartbasis

% Syntax:  v_rot = rodrigues(k,theta)

function v_rot = rodrigues_rot(k,theta)  
    k = k/sqrt(k(1)^2 + k(2)^2 + k(3)^2); % normalize rotation axis
    
    kk = 
      
end
