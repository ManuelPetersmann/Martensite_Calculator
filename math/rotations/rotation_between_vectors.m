function [R] = rotation_between_vectors(vec1,vec2,eps)

if nargin < 3
    eps = 1.e-4;
end
% call: rotation_between_vectors(vec1,vec2)
vec1 = vec1 / norm(vec1);
vec2 = vec2 / norm(vec2);
ang = acosd( dot ( vec1, vec2 ) );
if abs(ang) < eps
    %error('vectors are parallel')
    disp('Note: vectors are parallel - no rotation calculated');
    R = eye(3);
    return
end
ax = cross( vec1, vec2 ); % norming the vector happens in the function below
R = axis_angle_to_rotmat( ang, ax );

end

