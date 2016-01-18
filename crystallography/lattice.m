class BravaisLattice

function [ lattice_param ] = lattice_parameters()

bravais_type = input(prompt) % prompt all available lattices as strings

end

inPointGroup(R) - return true or false

inLatticeGroup(M) -"-


%    def isreversible(self)
%         
%         It is reversible if the symmetry group of ``U``
%         is a proper subgroup of the object's ``Laue`` group
% 
%         return whether the transformation is reversible
%         rtype boolean
%         raises AttributeError
% 
%                 if `U` has not been assigned
%         
%         if self.__U is None
%             raise AttributeError(U has not been initialized)
%         lg = Lattice(self.getU()).getLauegroup()
%         lg0 = self.getLaue()
%         return lg.order() = lg0.order() and lg0.hassubgroup(lg)        
%         if self.__U is None
%             raise AttributeError(U has not been initialized)
%         lg = Lattice(self.getU()).getLauegroup()
%         lg0 = self.getLaue()
%         return lg.order() = lg0.order() and lg0.hassubgroup(lg)


%     def __primitiveBase3D__(n, p):
%         '''
%         generate the base matrix for primitive cell of 3D Bravais lattices
%         '''
%         C = np.eye(3)
%         if n == 1:
%             # simple cubic
%             e1 = p[0] * np.array([1, 0, 0])
%             e2 = p[0] * np.array([0, 1, 0])
%             e3 = p[0] * np.array([0, 0, 1])
%         elif n == 2:
%             # face centered cubic
%             e1 = 0.5 * p[0] * np.array([1, 1, 0])
%             e2 = 0.5 * p[0] * np.array([0, 1, 1])
%             e3 = 0.5 * p[0] * np.array([1, 0, 1])
%             C = 0.5 * np.array([[1, 0, 1],
%                                 [1, 1, 0],
%                                 [0, 1, 1]])
%         elif n == 3:
%             # bcc
%             e1 = 0.5 * p[0] * np.array([1, 1, 1])
%             e2 = 0.5 * p[0] * np.array([-1, 1, 1])
%             e3 = 0.5 * p[0] * np.array([-1, -1, 1])
%             C = 0.5 * np.array([[1, -1, -1],
%                                 [1, 1, -1],
%                                 [1, 1, 1]])
%         elif n == 4:
%             # hexagonal
%             e1 = p[0] * np.array([1, 0, 0])
%             e2 = p[0] * np.array([0.5, np.sqrt(3) / 2, 0])
%             e3 = np.array([0, 0, p[1]])
%         elif n == 5:
%             # trigonal
%             #<111> is the 3fold axis
%             c = np.cos(p[1] * np.pi / 180)
%             a = p[0]
%             ty = np.sqrt((1 - c) / 6)
%             tz = np.sqrt((1 + 2 * c) / 3)
%             u = tz - 2 * np.sqrt(2) * ty
%             v = tz + np.sqrt(2) * ty
%             e1 = a / np.sqrt(3) * np.array([u, v, v])
%             e2 = a / np.sqrt(3) * np.array([v, u, v])
%             e3 = a / np.sqrt(3) * np.array([v, v, u])
%         elif n == 6:
%             # simple tetragonal
%             a = p[0]
%             c = p[1]
%             e1 = a * np.array([1, 0, 0])
%             e2 = a * np.array([0, 1, 0])
%             e3 = c * np.array([0, 0, 1])
%         elif n == 7:
%             # body centered tetragonal
%             a = p[0]
%             c = p[1]
%             e1 = (a / 2) * np.array([1, 1, c / a])
%             e2 = (a / 2) * np.array([-1, 1, c / a])
%             e3 = (a / 2) * np.array([-1, -1, c / a])
%             C = 0.5 * np.array([[1, -1, -1],
%                                 [1, 1, -1],
%                                 [1, 1, 1]])
%         elif n == 8:
%             # simple orthorhombic
%             a = p[0]
%             b = p[1]
%             c = p[2]
%             e1 = np.array([a, 0, 0])
%             e2 = np.array([0, b, 0])
%             e3 = np.array([0, 0, c])
%         elif n == 9:
%             # base centered orthorhombic
%             a = p[0]
%             b = p[1]
%             c = p[2]
%             e1 = np.array([a / 2, b / 2, 0])
%             e2 = np.array([-a / 2, b / 2, 0])
%             e3 = np.array([0, 0, c])
%             C = np.array([[0.5, -0.5, 0],
%                           [0.5, 0.5, 0],
%                           [0, 0, 1]])
%         elif n == 10:
%             # face centered orthorhombic
%             a = p[0]
%             b = p[1]
%             c = p[2]
%             e1 = np.array([a / 2, b / 2, 0])
%             e2 = np.array([0, b / 2, c / 2])
%             e3 = np.array([a / 2, 0, c / 2])
%             C = 0.5 * np.array([[1, 0, 1],
%                                 [1, 1, 0],
%                                 [0, 1, 1]])
%         elif n == 11:
%             # body centered orthorhombic
%             a = p[0]
%             b = p[1]
%             c = p[2]
%             e1 = np.array([a / 2, b / 2, c / 2])
%             e2 = np.array([-a / 2, b / 2, c / 2])
%             e3 = np.array([-a / 2, -b / 2, c / 2])
%             C = 0.5 * np.array([[1, -1, -1],
%                                 [1, 1, -1],
%                                 [1, 1, 1]])
%         elif n == 12:
%             # monoclinic unique axis b
%             a = p[0]
%             b = p[1]
%             c = p[2]
%             beta = radians(p[3])
%             e1 = np.array([a, 0, 0])
%             e2 = np.array([0, b, 0])
%             e3 = np.array([c * np.cos(beta), 0, c * np.sin(beta)])
%         elif n == 13:
%             # base centered monoclinic
%             a = p[0]
%             b = p[1]
%             c = p[2]
%             beta = radians(p[3])
%             e1 = np.array([a / 2, b / 2, 0])
%             e2 = np.array([-a / 2, b / 2, 0])
%             e3 = np.array([c * np.cos(beta), 0, c * np.sin(beta)])
%             C = np.array([[0.5, -0.5, 0],
%                           [0.5, 0.5, 0],
%                           [0, 0, 1]])
%         elif n == 14:
%             # triclinic
%             a = p[0]
%             b = p[1]
%             c = p[2]
%             alpha = radians(p[3])
%             beta = radians(p[4])
%             gamma = radians(p[5])
%             e1 = np.array([a, 0, 0])
%             e2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
%             e3 = np.array([c * np.cos(beta), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
%                            c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.cos(alpha) ** 2
%                                        - np.cos(beta) ** 2 - np.cos(gamma) ** 2) / np.sin(gamma)])
% 
%         return np.array([e1, e2, e3]).T, la.inv(C)

