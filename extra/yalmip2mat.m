function P=yalmip2mat(p,x)
% convert a yalmip polynomial to matrix
% Syntax: 
%   P=yalmip2mat(p,x);
%
% Inputs:
% p: yalmip polynomial
% x: sdpvar for polynomial variables
%
% Outputs:
% P.n: nb variables
% P.coef: list of coefficients
% P.pow: matrix of pows for monomials
%
% Author: Yi-Shuai NIU, 2019/08/13

[P.coef,monos]=coefficients(p,x);
k=length(P.coef); % number of monomials
P.n=length(x); % number of variables
P.pow=sparse(zeros(P.n,k));
for i=1:k
    P.pow(:,i)=degree(monos(i),x)';
end
end