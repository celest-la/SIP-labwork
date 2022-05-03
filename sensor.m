function [f,g,H,J] = sensor(x,params)
% x: value of the optimization variable
% params: structure containing the parameters of the objective function
% params.e = light intensity values
% params.r = sensor resistance values
% f,g,H,J: objective function, its gradient, hessian and Jacobian
%
E=params.e;
R=params.r/10^3;

f = sum(( R-x(1).*E.^(-x(2))).^2)/2;

g = [sum((R - x(1).*E.^(-x(2))).*(-E.^(-x(2)))) ; sum((R - x(1).*E.^(-x(2))).*(x(1).*E.^(-x(2)).*log(E))) ];

H(1,1)=sum(E.^(-2*x(2)));
H(1,2)=sum(R.*E.^(-x(2)).*log(E)-2*x(1).*E.^(-2*x(2)).*log(E));
H(2,1)=H(1,2);
H(2,2)=-sum(R.*x(1).*E.^(-x(2)).*((log(E)).^2))+sum(2.*x(1)^2.*E.^(-2*x(2)).*(log(E)).^2);


J(:,1)=-E'.^(-x(2));
J(:,2)=x(1).*E'.^(-x(2)).*log(E');