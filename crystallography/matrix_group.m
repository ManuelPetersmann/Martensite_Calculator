classdef matrix_group
    properties
        matrices
    end
    
    methods
        % constructor
        function obj = matrix_group( input )
            if isnumeric( input )
                obj.matrices = input;
            else
                obj.matrices = get_cubic_laue_group(); % TODO: write other cases e.g. hexagonal point group later
            end
        end
    end % end methods
end % end class

%-------------------------------------------------------------------
function cubic_laue_group = get_cubic_laue_group()
% cubic_point_group operations according to
% "Bhattacharya - Microstructure of Martensite" p.58
%
I = eye(3);
n = zeros(3,24);
alphas = zeros(1,24);
alphas(1) = 0; % angle and some vector to produce the identity first
n(:,1) = [1 0 0];
%
% 9 rotations around the x,y,z axes
counter = 1;
for i= 1:3
    alpha=90;
    for j = 1:3
        counter=counter+1;
        n(:,counter) = I(:,i);
        alphas(counter) = alpha;
        alpha=alpha+90;
    end
end
%
e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];
%
% the 6 -180 degree rotations around the face diagonals
alphas(11:16) = 180;
n(:,11) = 1/sqrt(2)*(e1+e2);
n(:,12) = 1/sqrt(2)*(e1-e2);
n(:,13) = 1/sqrt(2)*(e2+e3);
n(:,14) = 1/sqrt(2)*(e2-e3);
n(:,15) = 1/sqrt(2)*(e3+e1);
n(:,16) = 1/sqrt(2)*(e3-e1);
%
% the 8 rotations around the volume-diagonals
alphas(17:20) = 120;
n(:,17) = 1/sqrt(3)*(e1+e2+e3);
n(:,18) = 1/sqrt(3)*(e1+e2-e3);
n(:,19) = 1/sqrt(3)*(e1-e2+e3);
n(:,20) = 1/sqrt(3)*(-e1+e2+e3);
alphas(21:24) = 240;
n(:,21) = 1/sqrt(3)*(e1+e2+e3);
n(:,22) = 1/sqrt(3)*(e1+e2-e3);
n(:,23) = 1/sqrt(3)*(e1-e2+e3);
n(:,24) = 1/sqrt(3)*(-e1+e2+e3);
%
cubic_laue_group = zeros(3,3,24);
for i=1:size(n,2)
cubic_laue_group(:,:,i) = rot_originaxis_angle( alphas(i), n(:,i) );
end

end
%-------------------------------------------------------------------------
% TODO check for equivalence with above definitions (from PyStructTrans)
%
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
%-------------------------------------------------------------------------



