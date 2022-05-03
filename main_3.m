%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part of the code was used a first time to create the original and
% noised signals
%  x=simsignal(1000,'spectre');
%  x=x';

% figure(1)
% plot(x,'b')
% x0 = ones(length(x),1);
% params.n = length(x);
% N = params.n;
% 
% sigma = sqrt(sum(x(:).^2) * 10^(-RSB/10) / N);
% variance=sigma^2;
%  y = x + sigma*randn(size(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% This part of the code is used when x and y are saved %%%%%%%%%%%
load ex3_x.mat
load ex3_y.mat

x0 = ones(length(x),1);
params.n = length(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(22)
plot(x,'b')
figure(1)
plot(x,'b')

figure(2)
plot(y)

H=eye(length(x));
D=-eye(length(x)-1);
D=[D zeros(length(D),1)];
D=[zeros(length(D),1)'; D];
D=D+H;


params.y = y;
params.H = H;
params.D = D;


%%%%%%% enter all parameters and options %%%%%%%%%%%%%%%%%%%%%%%%%
params.TolX=10^-8;
params.TolF=10^-8;
params.MaxIter=10^3;
params.TolG=10^-8;
params.kmax=10000;
options.alpha=1;
options.method='Gauss-Newton'; 
%%Newton , Quasi-Newton, Gradient, Fletcher, Gauss-Newton, Levenberg
options.step = 'variable';
options.const=10^-3;
options.beta=0.75;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%descent method%%%%%%%%

%%calculation of x,f(lambda)
lambda=0:0.1:0.9;
error=zeros(length(lambda),1);
f0=zeros(length(lambda));
f1=zeros(length(lambda));
for i=1:length(lambda)
    params.lambda=lambda(i);
[xh]= optimdescent('func',params,options,x0);

%  figure(5)
%  subplot(length(lambda),1,i)
%  plot(xh)
[~,~,~,~,f0(i),f1(i)]=func(xh,params);
error(i)=norm(x-xh);
end

%%%determination of the minimal error
[mini,index]=min(error);
params.lambda=lambda(index);

[xh,result,xval]= optimdescent('func',params,options,x0);
figure(10)
plot(f1,f0)
xlabel ('f1(lambda)')
ylabel ('f0 (lambda)')
title ('evolution of f1(lambda) with the respect of f0(lambda)')
%%plot the best appromixation
figure (6)
 plot(xh)
 
 figure(22)
 hold on
 plot(xh,'g')
  hold off

%%Plot the error
figure (4)
plot(lambda, error)
title 'error'

%%Plot the gradient
figure (3)
semilogx(log10(result.grad))
title 'gradient norm'

%%Plot the f
figure (7)
semilogx(log10(result.crit))
title 'f'

figure(1)
params.lambda=0.2;
[xh]= optimdescent('func',params,options,x0);
hold on
 plot(xh,'k')
figure (1)
hold on
 params.lambda=1;
 [xh]= optimdescent('func',params,options,x0);
 plot(xh,'r')
 stralpha=sprintf('denoised(ideal alpha=%.1f)',lambda(index));
legend ('original signal', 'denoised (alpha=0.9)') 
hold off