function [ C_rot_full, C_rot_voig ] = rotateTensor4( R, C )
% call: rotateTensor4( R, C )
% rotates a fourth order tensor C_ijkl given a transformation matrix R
% in index notation this is written using Einsteins Summation convention as
% C_ijkl = R_im R_jn R_ko R_lp C_mnop 
% first output is the full 4-th order elastic tensor, second is the tensor
% in voigt notation (valid for fully symmetric 4th order tensor)

C_rot = zeros(3,3,3,3);

for i = 1:3
	for j = 1:3
		for k = 1:3
			for l = 1:3
				for m = 1:3
					for n = 1:3
						for o = 1:3
							for p = 1:3
								C_rot(i,j,k,l) = C_rot(i,j,k,l) + R(i,m) * R(j,n) * R(k,o) * R(l,p) * C(m,n,o,p);
                            end
                        end
                    end
                end
            end
        end
    end
end

C_rot_full = C_rot;

C_rot_voigt = [ C_rot(1,1,1,1) , C_rot(1,1,2,2) , C_rot(1,1,3,3) , C_rot(1,1,1,2) , C_rot(1,1,1,3) , C_rot(1,1,2,3) ;
                C_rot(2,2,1,1) , C_rot(2,2,2,2) , C_rot(2,2,3,3) , C_rot(2,2,1,2) , C_rot(2,2,1,3) , C_rot(2,2,2,3) ;
                C_rot(3,3,1,1) , C_rot(3,3,2,2) , C_rot(3,3,3,3) , C_rot(3,3,1,2) , C_rot(3,3,1,3) , C_rot(3,3,2,3) ;
                C_rot(1,2,1,1) , C_rot(1,2,2,2) , C_rot(1,2,3,3) , C_rot(1,2,1,2) , C_rot(1,2,1,3) , C_rot(1,2,2,3) ;
                C_rot(1,3,1,1) , C_rot(1,3,2,2) , C_rot(1,3,3,3) , C_rot(1,3,1,2) , C_rot(1,3,1,3) , C_rot(1,3,2,3) ;
                C_rot(2,3,1,1) , C_rot(2,3,2,2) , C_rot(2,3,3,3) , C_rot(2,3,1,2) , C_rot(2,3,1,3) , C_rot(2,3,2,3) ];

end