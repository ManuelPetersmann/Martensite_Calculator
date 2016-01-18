function [ cubic_point_group ] = get_cubic_pointgroup()

cubic_point_group = eye(3); % first enttrz is the identity

counter = 1;

    function[Rot] = rot(n, alpha)
alpha = alpha*(pi/180); 
n1=n(1);
n2=n(2);
n3=n(3);
         
Rot = round(  [(n1^2)*(1-cos(alpha))+cos(alpha)      n1*n2*(1-cos(alpha))-n3*sin(alpha)      n1*n3*(1-cos(alpha))+n2*sin(alpha)
              n2*n1*(1-cos(alpha))+n3*sin(alpha)    (n2^2)*(1-cos(alpha))+cos(alpha)         n2*n3*(1-cos(alpha))-n1*sin(alpha)
              n3*n1*(1-cos(alpha))-n2*sin(alpha)     n3*n2*(1-cos(alpha))+n1*sin(alpha)      (n3^2)*(1-cos(alpha))+cos(alpha)] ) ;
    end

% 9 rotations around the x,y,z axes
e= [1 0 0 
    0 1 0
    0 0 1];

for i= 1:3
    alpha=90;
    for j = 1:3
        counter=counter+1;
        cubic_point_group(:,:,counter) = rot(e(:,i), alpha);        
        alpha=alpha+90;
    end
end

% the 6 rotations around the face diagonals
e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];
alpha=180;

n = 1/sqrt(2)*(e1+e2);
cubic_point_group(:,:,11)= rotieren(n, alpha);

n = 1/sqrt(2)*(e1-e2);
cubic_point_group(:,:,12)= rotieren(n, alpha);

n = 1/sqrt(2)*(e2+e3);
cubic_point_group(:,:,13)= rotieren(n, alpha);

n = 1/sqrt(2)*(e2-e3);
cubic_point_group(:,:,14)= rotieren(n, alpha);

n = 1/sqrt(2)*(e3+e1);
cubic_point_group(:,:,15)= rotieren(n, alpha);

n = 1/sqrt(2)*(e3-e1);
cubic_point_group(:,:,16)= rotieren(n, alpha);

% the 8 rotations around the volume-diagonals
alpha = 120;

n = 1/sqrt(3)*(e1+e2+e3);
cubic_point_group(:,:,17)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1+e2-e3);
cubic_point_group(:,:,18)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1-e2+e3);
cubic_point_group(:,:,19)= rotieren(n, alpha);

n = 1/sqrt(3)*(-e1+e2+e3);
cubic_point_group(:,:,20)= rotieren(n, alpha);

alpha=240;

n = 1/sqrt(3)*(e1+e2+e3);
cubic_point_group(:,:,21)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1+e2-e3);
cubic_point_group(:,:,22)= rotieren(n, alpha);

n = 1/sqrt(3)*(e1-e2+e3);
cubic_point_group(:,:,23)= rotieren(n, alpha);

n = 1/sqrt(3)*(-e1+e2+e3);
cubic_point_group(:,:,24)= rotieren(n, alpha);

end

% CUBIC_LAUE_GROUP = MatrixGroup([
%     # identity
%     [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
%     # two fold rotations
%     [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # [100]
%     [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # [010]
%     [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],  # [001]
%     [[0, 1, 0], [1, 0, 0], [0, 0, -1]],  # [110]
%     [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],  # [1-10]
%     [[0, 0, 1], [0, -1, 0], [1, 0, 0]],  # [101]
%     [[0, 0, -1], [0, -1, 0], [-1, 0, 0]],  # [10-1]
%     [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],  # [011]
%     [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],  # [01-1]
%     # four fold rotations about [100]
%     [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
%     [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
%     # four fold rotations about [010]
%     [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
%     [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
%     # four fold rotations about [001]
%     [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
%     [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
%     # three-fold rotations about [111]
%     [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
%     [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
%     # three-fold rotations about [-111]
%     [[0, 0, -1], [-1, 0, 0], [0, 1, 0]],
%     [[0, -1, 0], [0, 0, 1], [-1, 0, 0]],
%     # three-fold rotations about [-1-11]
%     [[0, 1, 0], [0, 0, -1], [-1, 0, 0]],
%     [[0, 0, -1], [1, 0, 0], [0, -1, 0]],
%     # three-fold rotations about [1-11]
%     [[0, 0, 1], [-1, 0, 0], [0, -1, 0]],
%     [[0, -1, 0], [0, 0, -1], [1, 0, 0]]
% ])
% 
% __C1 = np.cos(np.pi / 3.0)
% __S1 = np.sin(np.pi / 3.0)
% HEX_LAUE_GROUP = MatrixGroup([
%     # identity
%     [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
%     # two fold rotations
%     [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # [100]
%     [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # [010]
%     [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],  # [001]
%     [[__C1, -__S1, 0], [__S1, __C1, 0],  [0, 0, 1]],
%     [[__C1, __S1, 0], [-__S1, __C1, 0],  [0, 0, 1]],
%     [[-__C1, -__S1, 0], [__S1, -__C1, 0],  [0, 0, 1]],
%     [[-__C1, __S1, 0], [-__S1, -__C1, 0],  [0, 0, 1]],
%     [[__C1, -__S1, 0], [-__S1, -__C1, 0],  [0, 0, -1]],
%     [[__C1, __S1, 0], [__S1, -__C1, 0],  [0, 0, -1]],
%     [[-__C1, -__S1, 0], [-__S1, __C1, 0],  [0, 0, -1]],
%     [[-__C1, __S1, 0], [__S1, __C1, 0],  [0, 0, -1]]
% ])


