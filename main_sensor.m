x1 = 0:0.01:4;
x2 = 0:0.01:2;
% 
params.r=[0.141 0.329 0.525 0.970 1.14 1.511 2.362 5.224 12.826 25.512];

params.e=[5013 2415 1558 1000 820 621 433 201 105 55];

N1 = length(x1);
N2 = length(x2);
F = zeros(N1,N2);

%matrix calculation
for i1 = 1:N1
    for i2 = 1:N2
        x = [x1(i1);x2(i2)];
        F(i1,i2) = sensor(x,params);
    end
end
figure(1)

meshc(x2,x1,log(F))


figure(2)

contour(x1,x2,log(F'),100)


%%%%%%% enter all parameters and options %%%%%%%%%%%%%%%%%%%%%%%%%
params.TolX=10^-8;
params.TolF=10^-8;
params.MaxIter=10^3;
params.TolG=10^-8;
params.kmax=10000;
params.n = 2;
options.alpha=1;
options.method='Newton'; 
%%Newton , Quasi-Newton, gradient, Fletcher, Gauss-Newton, Levenberg
options.step = 'variable';
options.const=10^-3;
options.beta=0.75;
x0=[1;1];
[ftest,gtest]=sensor(x0,params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%descent method%%%%%%%%

[xh,result,xval]= optimdescent('sensor',params,options,x0);

hold on

plot(xval(1,:),xval(2,:),'*r');

hold off

figure (3)
semilogx(log10(result.grad))

figure (4)
semilogx(log10(result.crit))