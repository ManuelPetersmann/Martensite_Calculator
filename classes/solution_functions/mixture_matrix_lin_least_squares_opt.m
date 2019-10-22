function x = mixture_matrix_lin_least_squares_opt( in )
% call: x = mixture_matrix_least_squares_opt( in )
% in... 3x3xn array of matrices
% optimize linear mixture rule x1 P1 + x2 P2 +.... xn Pn = Fcomp; 
% sum(xi)=1;  the Pi are the matrices in(:,:,i).
% https://de.mathworks.com/help/optim/ug/lsqlin.html
% Solves least-squares curve fitting problems of the form
% min || C*x - d|| with unequality, equality constraints and bounds!

nr = size(in,3);
% assemble C
Fc  = reshape(in,9,nr);
% assemble d
d = [1 0 0 0 1 0 0 0 1]';

%% constraints
Aeq = ones(1,nr); %[1 1];
beq = 1; % sum of fractions is one
% both parameters must be >= 0 !
A = (-1)*eye(nr);
%    [-1 0
%    0 -1];
b = zeros(1,nr); %[0 0]; % positivity constrained
lb = zeros(1,nr); %[0.;0.];
ub = ones(1,nr); %[1.;1.];

%% constrained linear least squares solution
%tic
x = lsqlin(Fc,d,A,b,Aeq,beq,lb,ub);
%toc

% H = Fc*x - d;
% H = reshape(H,3,3);
% % e.g. frobenius norm of displacement gradient
% d = frob(H);
% Fcomp = reshape(Fc*x,3,3);

%% quadratic programming solution
% min_x   1/2  x^T * H * x - f^t * x    with constraints
% see here why this linear least squares is also a quadratic programming
% problem (which can be solved more efficiently i guess...)
% https://math.stackexchange.com/questions/869204/are-constrained-linear-least-squares-and-quadratic-programming-the-same-thin

% H = Fc'*Fc;
% f = Fc' * d;
% tic
% x2 = quadprog(H,f,A,b,Aeq,beq,lb,ub)
% toc

end