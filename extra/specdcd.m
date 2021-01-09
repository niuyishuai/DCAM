function func = specdcd(x,Q,c)
% spectral dc decomposition to quadratic function
% f(x) = x'*Q*x + c'*x
% the function will return a dc function

[P,V]=eig(Q);
v=diag(V);
lstp=find(v>0);
lstn=find(v<0);

V1=V;
V1(lstn,lstn)=0;
Q1 = P*V1/P;

V2=V;
V2(lstp,lstp)=0;
Q2 = P*(-V2)/P;
if nargin == 3
    g= x'*Q1*x + c'*x;
    h= x'*Q2*x;
elseif nargin == 2
    g= x'*Q1*x;
    h= x'*Q2*x;
end

func = dcfunc(x,g,h);
end