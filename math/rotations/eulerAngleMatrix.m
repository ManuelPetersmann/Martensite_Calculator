% TODO specifiy other conventions, i.e. order of rotations and different -
% Roe, Kocks...
% axis using euler angles
function rot_mat = eulerAngleMatrix( alpha, beta, gamma, convention)
% This function takes the three euler Angles alpah, beta, gamma in degrees
% ( in order of rotation) and computes the rotation matrix with the
% convention of subsequent plane rotations  specified in the fourth argument.
% Available conventions are:
% 'bungee':  z, x', z'' 



if strcmpi(convention,'bungee')
    
    Rz = [cosd(gamma) sind(gamma) 0;
        -sind(gamma) cosd(gamma) 0;
        0         0      1];
    %
    Rx_strich = [ 1      0          0;
        0  cosd(beta) sind(beta);
        0 -sind(beta) cosd(beta)];
    %
    Rz_2strich = [ cosd(alpha)  sind(alpha)  0;
        -sind(alpha)  cosd(alpha)  0;
        0         0      1];
    %
    rot_mat = Rz * Rx_strich * Rz_2strich;
end

end
