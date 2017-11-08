function bool = in_O3( R )
% checks wheter a given matrix R is an orthogonal matrix
if (abs(R*R') - eye(3) < 1.e-5)
    if (abs( det(R) -1.) < 1.e-4)
        bool = true;
    else
        bool = false;
    end
else
    bool = false;
end
end
