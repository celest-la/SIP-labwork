function [f_lambda, g,Hessian,J,f0,f1] = func(x,params)

lambda = params.lambda;
y = params.y;
H = params.H;
D = params.D;

f0 = norm(y-H*x)^2;

f1 = norm(D*x)^2;

f_lambda = (1-lambda)*f0 + lambda*f1;

g = (1-lambda)*(-2*H'*(y-H*x)) + 2*lambda*(D'*D*x);

Hessian = 2*(1-lambda)*(H'*H)'+2*lambda*(D'*D)';

J = sqrt(2)*[-sqrt(1-lambda)*H;sqrt(lambda)*D];

end
