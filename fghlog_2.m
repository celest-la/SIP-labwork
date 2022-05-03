function[f,g,H,J]=fghlog(S,params)
beta=params.beta;
J=-sqrt(2)*params.K;

% if isempty(find(S<0,1))
f=(params.y-params.K*S)'*(params.y-params.K*S)+params.beta*(params.D*S)'*(params.D*S)-params.mu*norm(log(S))^2;
% else
% f=+Inf;
% end
g=-2*params.K'*(params.y-params.K*S)+2*beta*params.D'*(params.D*S)-params.mu*log(S)./S;
H=2*params.K'*params.K+2*beta*params.D'*params.D+2*params.mu*diag(1-log(S)./(S.^2));
end
