% TODO specifiy other conventions, i.e. order of rotations and different
% axis using euler angles
function rot_mat = eulerAngleMatrix( alpha, beta, gamma, convention)
% This function takes the three euler Angles alpah, beta, gamma  
% ( in order of rotation) and computes the rotation matrix with the
% convention of subsequent plane rotations  specified in the fourth argument.
% Available conventions are:
% 'bungee':  z, x', z'' 


if strcmpi(convention,'bungee')
    
    Rz = [cos(gamma) sin(gamma) 0;
        -sin(gamma) cos(gamma) 0;
        0         0      1];
    %
    Rx_strich = [ 1      0          0;
        0  cos(beta) sin(beta);
        0 -sin(beta) cos(beta)];
    %
    Rz_2strich = [ cos(alpha)  sin(alpha)  0;
        -sin(alpha)  cos(alpha)  0;
        0         0      1];
    %
    rot_mat = Rz * Rx_strich * Rz_2strich;
end

end
