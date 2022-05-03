function [f,g,H,J] = rosenbrock(x,params)
% x : value of the optimization variable
% f,g,h : objective function value, gradient and Hessian
% params : structure containing the objective function parameters
% params.n : value of n
% params.b : value of b
%
n=params.n;
r=zeros(2*(n-1),1);
g=zeros(n,1);
b=params.b;
for i=1:n-1
    r(i)=sqrt(params.b)*(x(i+1)-x(i)^2);
    r(i*2)=1-x(i);
end
   f=sum(r.^2); 
  
   %%calculation of gradient
    g(1) = -4*params.b*x(1)*(x(2)-x(1)^2)-2*(1-x(1));
    g(n)=2*params.b*(x(n)-x(n-1)^2);
    
  if n>2
    for i=2:n-1
      g(i)=2*b*((x(i)^2-x(i-1))*2*x(i)+(x(i)-x(i-1)^2))-2*(1-x(i));
    end
    end


%%calcuation of the Hessian
H(1,1)=b*(12*x(1)^3-4*x(2))+2;
H(n,n)=2*b;
H(n,n-1)=-4*b*x(n-1);
H(n-1,n)=-4*b*x(n-1);
  
  if n>2
    for i=2:n-1
      H(i,i-1)=-4*b*x(i-1);
      H(i-1,i)=-4*b*x(i-1);
      H(i,i)=-4*b*x(i+1)+12*b*x(i)^2*x(i)+2;
    end
  end

  %%calculation of the Jacobian
  J=zeros(2*(n-1),2);
  for i=1:n-1
    J(i,1)=-2*sqrt(params.b)*x(i);
    J(i,2)=sqrt(params.b);
    J(2*i,1)=-1;
    
  end

end
 