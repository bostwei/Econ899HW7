function [w] = utility(inputPar,kk)
%UTILITY calculate the utiluty function given k and K
%   Detailed explanation goes here

global alpha beta 

k = inputPar.k;
K = inputPar.K;
A = inputPar.A;
B = inputPar.B;

% calculate the r and w 
w = (1-alpha) * K.^alpha;
r = alpha * K.^(alpha-1);

% calculate the consumption
c = r.*k + w - kk';
c(c<0) = 0;
u = log(c);

v = A + B * log(kk);

w = u + beta * v';
end

