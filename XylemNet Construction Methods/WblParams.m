function [A,B] = WblParams(mu,sig)

f = @(k)sig^2/mu^2-gamma(1+2./k)./gamma(1+1./k).^2+1;
B = fzero(f,[0.1 200]);       % Solve for shape parameter
A = mu/gamma(1+1/B); % Substitue to find scale parameter

if abs(B-0.1)<0.01 || abs(B-200)<0.01
   warning('Weibull parameters might not have been found') 
end
end