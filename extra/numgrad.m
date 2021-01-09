function fgrad=numgrad(F,x0,delta)
% compute numerical gradient for any shape of x0 using formulation
% dF/dx (x0) = (F(x+delta) - F(x-delta))/(2*delta)
% fgrad=numgrad(F,x0)
% F: a function handle
% x0: a given point
% delta: step size
if nargin < 3
    delta=0.01;
end
n=numel(x0);
fgrad=zeros(size(x0));
for i=1:n
    x=x0;
    % compute fx1=f(x0+delta*ei)
    x(i)=x(i)+delta;
    fx1=F(x);
    % compute fx2=f(x0-delta*ei)
    x(i)=x(i)-2*delta;
    fx2=F(x);
    % compute (fx1-fx2)/(2*delta)
    fgrad(i)=(fx1-fx2)/(2*delta);
end
end