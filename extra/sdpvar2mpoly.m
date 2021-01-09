function mpolyf = sdpvar2mpoly(yalmipf,yalmipx)
n=length(yalmipx); % number of variables
mpolyf=mpolymat(n,size(yalmipf));
for idx=1:numel(yalmipf)
    % convert yalmip object to mpoly object
    [coef,monos]=coefficients(yalmipf(idx),yalmipx);
    k=length(coef); % number of monomials
    mpolyf(idx).k=k;
    mpolyf(idx).coef=coef; % k*1 vector
    mpolyf(idx).pow=zeros(k,n); % k*n matrix
    for i=1:k
        mpolyf(idx).pow(i,:)=degree(monos(i),yalmipx);
    end
end
end