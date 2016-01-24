% This file provides several functions for GEOMETRIC TRANSFORMATIONS 
% such as rotating (active rotation), mirroring, scaling, shearing

function funcs = geom_transforms( input_args )

% Idea -> use this function to calculate correspondance matrix based on 
% description in "Bhadeshia - Worked examples in the geometry of crystals"
Gauss =  linalg::gaussElim(A)


end

%TODO rewrite to Matlab
%  def to_primitive_Miller(self, idx):
%         """
%         convert a Miller index or a list of Miller indices
%         from conventional base to primitive base
% 
%         :param idx: a Miller index or a list of Miller indices
%         :return: the Miller index or the list of Miller indices in primitive base
%         """
%         try:
%             idx = np.array(idx)
%             if idx.ndim == 1:
%                 idx = np.array([idx])
%             elif idx.ndim > 2:
%                 raise ValueError("invalid input")
%             if idx.shape[1] is not self.dimension():
%                 raise ValueError("dimensions of lattice and indices not match")
%         except:
%             raise ValueError("invalid input of index (indices)")
% 
%         res = (self.__C.dot(idx.T)).T
%         return res[0] if len(idx) == 1 else res
% 
%     def toPrimitive(self, idx):
%         warn("toPrimitive() is deprecated, please use to_primitive() instead.", DeprecationWarning)
%         return self.to_primitive()
% 
%     def to_conventional_Miller(self, idx):
%         """
%         convert a Miller index or a list of Miller indices
%         from primitive to conventional
% 
%         :param idx: a Miller index or a list of Miller indices
%         :return: the Miller index or the Miller list of indices in conventional base
%         """
%         try:
%             if not isinstance(idx, np.ndarray):
%                 idx = np.array(idx)
%             if idx.ndim == 1:
%                 idx = np.array([idx])
%             elif idx.ndim > 2:
%                 raise ValueError("invalid input")
%             if idx.shape[1] is not self.dimension():
%                 raise ValueError("dimensions of lattice and indices not match")
%         except:
%             raise ValueError("invalid input of index (indices)")
% 
%         res = (la.inv(self.__C).dot(idx.T)).T
%         return res[0] if len(idx) == 1 else res
