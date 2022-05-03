function [f,g, H, J] = func2(S,params)

y = params.y;
K = params.K;

mu = params.mu;
lambda = params.lambda;
D = params.D;

W = eye(length(S));
    for i= 1:length(S)
        if(S(i)>=0)
           W(i,i) = 0; 
        end
    end

f = norm(y - K*S)^2 + lambda*norm(D*S)^2 + mu*norm(W*S)^2;
g = (-2*K'*(y-K*S)) + 2*lambda*(D'*D*S) +  2*mu*(W'*W*S);
H = 2*(K'*K)'+2*lambda*(D'*D)' + 2*mu*(W'*W)';
J = sqrt(2)*[-K;sqrt(lambda)*D];
end