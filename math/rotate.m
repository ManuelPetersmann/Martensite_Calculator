function out = rot(input, rot_mat)
if size(input) == [1, 3] %input must be a colon vector
    % Note that generally instead of the transpose here is an inverse, but
    % since for a rotatoin R*R'=I holds its the same
    out = rot_mat' * input; % passive rotation
elseif size(input) == [3, 3]
    out = zeros(3);
    for i = 1:3
        out(:,i) = rot_mat' * input(:,i);
    end
elseif size(input) == [3, 3, 3, 3]
    % TODO implement for 4th order tensor
end
end


