%initialization
x1 = -4:0.05:4;
x2 = -5:0.05:20;
params.b = 2;
params.n = 2;
N1 = length(x1);
N2 = length(x2);
F = zeros(N1,N2);

%matrix calculation
for i1 = 1:N1
    for i2 = 1:N2
        x = [x1(i1);x2(i2)];
        F(i1,i2) = rosenbrock(x,params);
    end
end
figure(1)
meshc(x2,x1,F)
title '3D visualization of the Rosenbrock function'

figure (2)
contour(x2,x1,F,100)
hold on
title 'visualisation of the Rosenbrock function with contour'

%%%%%%% enter all parameters and options %%%%%%%%%%%%%%%%%%%%%%%%%
params.TolX=10^-8;
params.TolF=10^-8;
params.MaxIter=10^3;
params.TolG=10^-8;
params.kmax=10000;
options.alpha=1;
options.method='gradient';
%%gradient, Newton , Quasi-Newton, Fletcher, Gauss-Newton, Levenberg
options.step = 'variable';
options.const=10^-3;
options.beta=0.75;
x0=[-2;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%descent method%%%%%%%%

[xh,result,xval]= optimdescent('rosenbrock',params,options,x0);

plot(xval(2,:),xval(1,:),'*');

hold off
figure (3)
semilogx(log10(result.grad),'b')
legend 'Gradient descent'
title 'Gradient=f(iteration)'

figure (4)
semilogx(log10(result.crit),'b')

title 'Criterion=f(iteration)'
time.gradient=result.time;
options.method='Newton';
[xh,result,xval]= optimdescent('rosenbrock',params,options,x0);

figure (3)
hold on
semilogx(log10(result.grad),'r')
hold off

figure (4)
hold on
semilogx(log10(result.crit),'r')
hold off
time.newton=result.time;
options.method='Quasi-Newton';
[xh,result,xval]= optimdescent('rosenbrock',params,options,x0);

figure (3)
hold on
semilogx(log10(result.grad),'g')
hold off

figure (4)
hold on
semilogx(log10(result.crit),'g')
hold off
time.quasinewton=result.time;
options.method='Fletcher';
[xh,result,xval]= optimdescent('rosenbrock',params,options,x0);

figure (3)
hold on
semilogx(log10(result.grad),'y')
hold off

figure (4)
hold on
semilogx(log10(result.crit),'y')
hold off
time.fletcher=result.time;

options.method='Gauss-Newton';
[xh,result,xval]= optimdescent('rosenbrock',params,options,x0);

figure (3)
hold on
semilogx(log10(result.grad),'c')
hold off

figure (4)
hold on
semilogx(log10(result.crit),'c')
hold off
time.gauss=result.time;
options.method='Levenberg';
[xh,result,xval]= optimdescent('rosenbrock',params,options,x0);

figure (3)
hold on
semilogx(log10(result.grad),'m')

legend ('Gradient Descent', 'Newton', 'Quasi-Newton', 'Fletcher-Reeves','Gauss-Newton', 'Levenberg')
hold off

figure (4)
hold on
semilogx(log10(result.crit),'m')

legend ('Gradient Descent', 'Newton', 'Quasi-Newton', 'Fletcher-Reeves','Gauss-Newton', 'Levenberg')
hold off
time.levenberg=result.time;
figure(3)
