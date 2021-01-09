function diagnostics=callconvexsolver_bis(f, l, C, x0, solver, options)
% f: the convex part of objective function
% l: the linear coefs of objective function
% C: constraints
% x0: initial point
% solver: name of solver, fmincon|ipopt|knitro|filterSD|wolfe
% options: option for solver
if isempty(l)
    l=zeros(size(x0));
end
if nargin < 6
    options = optimset('Algorithm', 'interior-point', 'Display','off');
end

if strcmpi(solver,'gd')
    funcs.objective = @(x)f.eval(x) + l'*x;
    funcs.gradient = @(x)numgrad(funcs.objective,x);
    alpha = 0.1;
    delta = inf;
    while (delta > 1e-3)
        x1=BPPPA(funcs.gradient(x0),1);
        d = alpha*(x1-x0);
        x0 = x0 + d;
        delta = norm(d);
        if alpha >1e-3
            alpha = alpha/2;
        end
    end
    diagnostics.x=x0;
    diagnostics.fval=funcs.objective(x0);
    diagnostics.problem=0;
    return;
end


if strcmpi(solver,'pg-bpppa')
    funcs.objective = @(x)f.eval(x) + l'*x;
    funcs.gradient = @(x)numgrad(funcs.objective,x);
    delta = inf;
    while (delta > 1e-7)
        x1=BPPPA(funcs.gradient(x0),1);
        delta = norm(x1-x0);
        x0=x1;
    end
    diagnostics.x=x1;
    diagnostics.fval=funcs.objective(x1);
    diagnostics.problem=0;
    return;
end

if strcmpi(solver,'wolfe')
    % funcs structure
    funcs.objective = @(x)f.eval(x) + l'*x;
    funcs.gradient = @(x)numgrad(funcs.objective,x);
    %funcs.gradient = @(x)df.eval(x)+l;
    if checkfeas(x0)==false
        x0 = knitromatlab(funcs.objective,x0,[],[],C.Aeq,C.beq,C.lb,C.ub,[],[],options);
        diagnostics.x=x0;
        diagnostics.fval=funcs.objective(x0);
        diagnostics.problem=0;
        return;
    end
    n=length(l);
    ub=ones(n,1);
    while true
        grad=funcs.gradient(x0);
        lb=-ub;
        for i=1:n
            if abs(x0(i))<=1e-6
                lb(i)=0;
            end
        end
        d=cplexlp(grad,[],[],C.Aeq,0,lb,ub);
        if norm(d)<1e-6
            diagnostics.x=x0;
            diagnostics.fval=funcs.objective(x0);
            diagnostics.problem=0;
            return;
        end
        % line search in descent direction
        linfunc = @(s)funcs.objective(x0+s*d);
        idx=find(d<0);
        if isempty(idx)
            alpha=inf;
        else
            alpha = min(-x0(idx)./d(idx));
        end
        s=fminbnd(linfunc,0,alpha);
        f0=funcs.objective(x0);
        x0=x0+s*d;
        f1=funcs.objective(x0);
        if abs(f1-f0) < 1e-6
            diagnostics.x=x0;
            diagnostics.fval=funcs.objective(x0);
            diagnostics.problem=0;
            return;
        end
        
        %fprintf('the optimal value is : %.8e\n',funcs.objective(x0));
    end
end

if strcmpi(solver,'wolfe1')
    % funcs structure
    funcs.objective = @(x)f.eval(x) + l'*x;
    funcs.gradient = @(x)numgrad(funcs.objective,x);
    %funcs.gradient = @(x)df.eval(x)+l;
    grad=funcs.gradient(x0);
    x=cplexlp(grad,[],[],C.Aeq,C.beq,C.lb,C.ub,x0);
    if checkfeas(x0)==false % 初始点不可行
        diagnostics.x=x;
        diagnostics.fval=funcs.objective(x);
        diagnostics.problem=0;
        return;
    end
    d=x-x0; % descent direction
    if grad'*d >= 1e-6
        diagnostics.x=x0;
        diagnostics.fval=funcs.objective(x0);
        diagnostics.problem=0;
        return;
    end
    % line search in descent direction
    f0 = funcs.objective(x0);
    f1 = funcs.objective(x);
    if f1 < f0
        diagnostics.x=x;
        diagnostics.fval=f1;
        diagnostics.problem=0;
        return;
    end
    linfunc = @(s)funcs.objective(x0+s*d);
    [s,fval]=fminbnd(linfunc,0,1);
    x=x0+s*d;
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=0;
    return;
end

if strcmpi(solver,'wolfe2')
    % funcs structure
    funcs.objective = @(x)f.eval(x) + l'*x;
    funcs.gradient = @(x)numgrad(funcs.objective,x);
    %funcs.gradient = @(x)df.eval(x)+l;
    grad=funcs.gradient(x0);
    x=cplexlp(grad,[],[],C.Aeq,C.beq,C.lb,C.ub,x0);
    if checkfeas(x0)==false % 初始点不可行
        diagnostics.x=x;
        diagnostics.fval=funcs.objective(x);
        diagnostics.problem=0;
        return;
    end
    d=x-x0; % descent direction
    if grad'*d >= 1e-6
        diagnostics.x=x0;
        diagnostics.fval=funcs.objective(x0);
        diagnostics.problem=0;
        return;
    end
    % line search in descent direction
    f0 = funcs.objective(x0);
    beta = 0.618;
    alpha = 1;
    while (alpha>1e-5)
        x=x0+alpha*d;
        f1=funcs.objective(x);
        if f1 < f0 % x is always feasible since [x,x0] in C
            diagnostics.x=x;
            diagnostics.fval=f1;
            diagnostics.problem=0;
            return;
        end
        alpha=beta*alpha;
    end
    diagnostics.x=x0;
    diagnostics.fval=f0;
    diagnostics.problem=0;
    return;
end

if strcmpi(solver,'knitro')
    func = @(x)f.eval(x) + l'*x;
    [x,fval,exitflag] = ...
        knitromatlab(func,x0,[],[],C.Aeq,C.beq,C.lb,C.ub,[],[],options);
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=exitflag;
    return;
end

if strcmp(solver,'fmincon')
    func = @(x)f.eval(x) + l'*x;
    [x,fval,exitflag] = ...
        fmincon(func,x0,[],[],C.Aeq,C.beq,C.lb,C.ub,[],options);
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=exitflag-1;
    return;
end

if strcmpi(solver,'ipopt')
    % funcs structure
    funcs.objective = @(x)f.eval(x) + l'*x;
    %grad=f.jacobian' + l;
    %funcs.gradient = @(x)grad.eval(x);
    funcs.gradient = @(x)numgrad(funcs.objective,x);
    funcs.constraints = @(x)C.Aeq*x;
    jac=sparse(ones(1,f.n));
    funcs.jacobian = @(x)jac;
    funcs.jacobianstructure = funcs.jacobian;
    % opts structure
    opts.lb = C.lb;
    opts.ub = C.ub;
    opts.cl = C.beq;
    opts.cu = C.beq;
    opts.ipopt = ipoptset;
    opts.ipopt.hessian_approximation = 'limited-memory';
    opts.ipopt.mu_strategy           = 'adaptive';
    opts.ipopt.tol                   = 1e-7;
    opts.ipopt.print_level           = 0;
    
    [x,output] = ...
        ipopt(x0,funcs,opts);
    diagnostics.x=x;
    diagnostics.fval=f.eval(x);
    diagnostics.problem=output.status;
    return;
end

if strcmpi(solver,'filterSD')
    % functions
    func = @(x)f.eval(x) + l'*x;
    grad = @(x)numgrad(func,x);
    nlcon = @(x)C.Aeq*x;
    nljac = @(x)ones(1,f.n);
    % opts structure
    opts = optiset('solver','filterSD');
    opts.tolrfun = 1e-7;
    opts.tolafun = 1e-7;
    
    [x,fval,exitflag] = ...
        opti_filtersd(func, grad, x0, C.lb, C.ub, nlcon, nljac, [], C.beq, C.beq, opts);
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=exitflag-1;
    return;
end

if strcmp(solver,'baron')
    % function handle (need to build function from string, can not use f.eval)
    str=f.func2str;
    str=strcat('@(x)',str);
    func = str2func(str);
    A = C.Aeq;
    rl = C.beq;
    ru = C.beq;
    lb = zeros(f.n,1);
    ub = ones(f.n,1);
    % opts structure
    opts = baronset;
    opts.prlevel = 0;
    opts.EpsR = 1e-9;
    opts.EpsA = 1e-6;
    
    [x,fval,exitflag] = ...
        baron(func,A,rl,ru,lb,ub,[],[],[],[],[],opts);
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=exitflag-1;
    return;
end

error('no solver found');

end

function r=checkfeas(x)
%min(check(obj.dcp.C))>=-1e-7
if abs(sum(x)-1) <=1e-7 && min(x) >= -1e-7
    r=true;
else
    r=false;
end
end

