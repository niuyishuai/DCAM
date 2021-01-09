function p=mat2yalmip(P,x)
% convert a matrix format polynomial to yalmip object
% Syntax: 
%   p=mat2yalmip(P,x);
%
% Inputs:
% P: matrix format polynomial structure
%  -P.n: nb variables
%  -P.coef: list of coefficients
%  -P.pow: matrix of pows for monomials
% x: sdpvar for polynomial variables
%
% Outputs:
% p: yalmip polynomial
%
% Author: Yi-Shuai NIU, 2019/08/13

if numel(x)~=P.n
    error('dimensions of P and x do not match!');
end
k=size(P.pow,2); % number of monomials
p=prod((x(:)*ones(1,k)).^P.pow)*P.coef(:);
end