function[f,g,H,J]=fghlog(S,params)
beta=params.beta;
J=-sqrt(2)*params.K;

% if isempty(find(S<0,1))
% f=(params.y-params.K*S)'*(params.y-params.K*S)+params.beta*(params.D*S)'*(params.D*S)-params.mu*sum(log(S));
% else
% f=+Inf;
% end
% g=-2*params.K'*(params.y-params.K*S)+2*beta*params.D'*(params.D*S)-params.mu./S;
% H=2*params.K'*params.K+2*beta*params.D'*params.D+diag(params.mu./(S.^2));
if isempty(find(S<0,1))
logS=log(S);
else
    logS=+Inf;
end
%if isempty(find(S<0,1))
f=(params.y-params.K*S)'*(params.y-params.K*S)+params.beta*(params.D*S)'*(params.D*S)+params.mu*norm(logS)^2;
% else
%  f=+Inf;
% end
g=-2*params.K'*(params.y-params.K*S)+2*beta*params.D'*(params.D*S)-params.mu*logS./S;
H=2*params.K'*params.K+2*beta*params.D'*params.D+2*params.mu*diag(1-logS./(S.^2));
end

