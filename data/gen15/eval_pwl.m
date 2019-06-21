function [x,y] = eval_pwl(q,c,qmin,qmax,nPts)

% assert(abs(qmin) > 0.01); % ensure generator min level is positive

% x = linspace(qmin,qmax,nPts-1)';
x = linspace(qmin,qmax,nPts)';

[~,k] = histc(x,q); % vectorized binning
tol = 0.0001
k(abs(k - length(q)) < tol) = length(q) - 1;
t = (x - q(k))./(q(k+1) - q(k));
y = (1-t).*c(k) + t.*c(k+1);

% return with concatenating 0 to the power levels
% x = [0;x];
% y = [0;y];
