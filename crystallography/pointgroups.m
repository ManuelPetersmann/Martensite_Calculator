function funcs = pointgroups()
funcs.determine_pointgroup = @determine_pointgroup;
%funcs.inPointGroup = @in_pointgroup;
end


function pointgroup = determine_pointgroup( Bravais_type )
% For martensite transformations the point_group of the parent phase is
% important (mostly cubic or hexagonal, implement others if needed...)
if strcmpi(Bravais_type, 'cubic')
    pointgroup = get_cubic_pointgroup()
end
if strcmpi(Bravais_type, 'hexagonal')
% TODO implement        pointgroup = get_cubic_pointgroup()
end
end

    
function [ cubic_point_group ] = get_cubic_pointgroup()
% cubic_point_group operations according to 
% "Bhattacharya - Microstructure of Martensite" p.58
%
r = coord_transforms;
%
cubic_point_group = eye(3); % first entry is the identity
%
% 9 rotations around the x,y,z axes
e= [1 0 0 
    0 1 0
    0 0 1];
%
counter = 1;
for i= 1:3
    alpha=90;
    for j = 1:3
        counter=counter+1;
        cubic_point_group(:,:,counter) = r.rot_mat_axis(e(:,i), alpha);        
        alpha=alpha+90;
    end
end
%
% the 6 rotations around the face diagonals
e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];
alpha=180;
%
n = 1/sqrt(2)*(e1+e2);
cubic_point_group(:,:,11)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(2)*(e1-e2);
cubic_point_group(:,:,12)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(2)*(e2+e3);
cubic_point_group(:,:,13)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(2)*(e2-e3);
cubic_point_group(:,:,14)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(2)*(e3+e1);
cubic_point_group(:,:,15)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(2)*(e3-e1);
cubic_point_group(:,:,16)= r.rot_mat_axis(n, alpha);

% the 8 rotations around the volume-diagonals
alpha = 120;

n = 1/sqrt(3)*(e1+e2+e3);
cubic_point_group(:,:,17)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(3)*(e1+e2-e3);
cubic_point_group(:,:,18)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(3)*(e1-e2+e3);
cubic_point_group(:,:,19)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(3)*(-e1+e2+e3);
cubic_point_group(:,:,20)= r.rot_mat_axis(n, alpha);

alpha=240;

n = 1/sqrt(3)*(e1+e2+e3);
cubic_point_group(:,:,21)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(3)*(e1+e2-e3);
cubic_point_group(:,:,22)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(3)*(e1-e2+e3);
cubic_point_group(:,:,23)= r.rot_mat_axis(n, alpha);

n = 1/sqrt(3)*(-e1+e2+e3);
cubic_point_group(:,:,24)= r.rot_mat_axis(n, alpha);
end

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
         %--------------------------------------------------------------
%         TODO - rewrite to matlab
%         function[bool] = inPointGroup(Q) 
%         Check if `Q` is in the point group of the lattice.
%         This method works for lattices in any dimension.
%         if Q is not an orthogonal matrix, return false
%         try:
%             Q = np.array(Q)
%         except Exception:
%             return False
%         if not _in_O3(Q):
%             return False
%         if len(Q) != self.__N:
%             return False
%         return self == Lattice(np.dot(Q, self.__E))
         %--------------------------------------------------------------
         % TODO - rewrite to matlab
%         function[bool] = in lattice_group --> find out what they mean by
%         lattice group, maybe e.g. tetragonal containing P,I...
%         """
%         Check if `M` is in the lattice group of the lattice.
%         This method works for lattices in any dimension.
%         """
%         # if Q is not an orthogonal matrix, return false
%         try:
%             M = np.array(M)
%         except Exception:
%             return False
%         M_int = np.rint(M)
%         if not np.allclose(M, M_int, atol=1E-8):
%             return False
%         # QE = EM
%         E = self.base()
%         Q = E.dot(M).dot(la.inv(E))
%         return _in_O3(Q)

