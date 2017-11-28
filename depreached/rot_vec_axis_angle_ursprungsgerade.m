function[ v_out, Rot ] = rot_vec_axis_angle_ursprungsgerade(v_in, alpha, n )
% call as: rot_vec_axis_angle(v_in, alpha, n )
% v_in - (row) vector to rotate
% alpha - angle to rotate
% n - around which should be rotated (right-handed system, positive sense)
alpha = deg2rad(alpha); %umrechnen von Grad aus Argument auf Radianten für cos und sin
n1=n(1) / norm(n,2);
n2=n(2) / norm(n,2);
n3=n(3) / norm(n,2);

Rot = [(n1^2)*(1-cos(alpha))+cos(alpha)      n1*n2*(1-cos(alpha))-n3*sin(alpha)      n1*n3*(1-cos(alpha))+n2*sin(alpha)
    n2*n1*(1-cos(alpha))+n3*sin(alpha)   (n2^2)*(1-cos(alpha))+cos(alpha)         n2*n3*(1-cos(alpha))-n1*sin(alpha)
    n3*n1*(1-cos(alpha))-n2*sin(alpha)   n3*n2*(1-cos(alpha))+n1*sin(alpha)       (n3^2)*(1-cos(alpha))+cos(alpha)];

v_out = Rot * v_in';

end