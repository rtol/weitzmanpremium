function [a LL] = fitploeg(temp,imp)
%function [a LL] = fitploeg(temp,imp)
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.1-matlab-global
%
%This function is part of FUND 4.1 MG
%It fits Van der Ploeg's impact function
%
%Richard Tol, 18 March 2018
%This code is protected by the MIT License

X = exp(temp)-1;
a = inv(X'*X)*X'*imp;
res = imp - X*a;
ssr = res'*res;
n = length(imp);
p = 1;
ssq = ssr/(n-p);

LL= -0.5*n*log(2*pi*ssq) - 0.5*ssr/ssq;