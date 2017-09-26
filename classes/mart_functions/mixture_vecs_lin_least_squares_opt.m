function [x,eps_net] = mixture_vecs_lin_least_squares_opt( in )
% call: [x,d] = mixture_vecs_lin_least_squares_opt( in )
% in... 3xn array of colum vectors
% optimize linear mixture rule x1 vec1 + x2 vec2 +.... xn vecn; sum(x_i)=1
% https://de.mathworks.com/help/optim/ug/lsqlin.html
% Solves least-squares curve fitting problems of the form
% min || C*x - d|| with unequality, equality constraints and bounds!

nr = size(in,2);
% assemble d
d = [0 0 0]';

%% constraints
Aeq = ones(1,nr); %[1 1];
beq = 1; % sum of fractions is one
% both parameters must be >= 0 !
A = (-1)*eye(nr);
%    [-1 0
%    0 -1];
b = zeros(1,nr); % positivity constrained
lb = zeros(1,nr); 
ub = ones(1,nr); 

%% constrained linear least squares solution
%tic
x = lsqlin(in,d,A,b,Aeq,beq,lb,ub);
%toc

eps_net = norm( in*x );

end