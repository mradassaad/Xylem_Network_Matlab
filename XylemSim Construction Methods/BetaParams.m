function [A,B] = BetaParams(mu,sig)

t = mu/(1-mu);
f = @(B)sig^2 - t*B^2/(B^2*(1-t)^2*(B*(1+t)+1));
B = fzero(f,[0.1 200]); % Solve for shape parameter B
A = t*B; % Substitue to find shape parameter A

end