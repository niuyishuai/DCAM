function P=sym2mat(p,x)
% convert a symbolic polynomial to matrix
% Syntax: 
%   P=sym2mat(p,x);
%
% Inputs:
% p: symbolic polynomial
% x: symbolic polynomial variables
%
% Outputs:
% P.n: nb variables
% P.coef: list of coefficients
% P.pow: matrix of pows for monomials
%
% Author: Yi-Shuai NIU, 2019/08/19

[P.coef,monos]=coeffs(p,x);
k=length(P.coef); % number of monomials
P.n=length(x); % number of variables
P.pow=sparse(zeros(P.n,k));
for i=1:k
    for j=1:P.n
        P.pow(j,i)=polynomialDegree(monos(i),x(j));
    end
end
end