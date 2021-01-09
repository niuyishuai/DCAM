function diagnostics=callconvexsolver(f, C, x0, solver, options)
% f: objective function
% C: constraints
% x0: initial point
% solver: name of solver, fmincon|ipopt|gurobi|cplex|knitro|filterSD|CVX
% options: option for solver
if nargin < 5
    options = optimset('Algorithm', 'interior-point', 'Display','off');
end
if strcmp(solver,'gurobi')
    func = @(x)f.eval(x);
    gurobi
    [x,fval,exitflag] = ...
        knitromatlab(func,x0,[],[],C.Aeq,C.beq,C.lb,C.ub,[],[],options);
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=exitflag;
    return;
end
if strcmp(solver,'knitro')
    func = @(x)f.eval(x);
    [x,fval,exitflag] = ...
        knitromatlab(func,x0,[],[],C.Aeq,C.beq,C.lb,C.ub,[],[],options);
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=exitflag;
    return;
end

if strcmp(solver,'fmincon')
    func = @(x)f.eval(x);
    [x,fval,exitflag] = ...
        fmincon(func,x0,[],[],C.Aeq,C.beq,C.lb,C.ub,[],options);
    diagnostics.x=x;
    diagnostics.fval=fval;
    diagnostics.problem=exitflag-1;
    return;
end

if strcmp(solver,'ipopt')
    % funcs structure
    funcs.objective = @(x)f.eval(x);
    grad=f.jacobian';
    funcs.gradient = @(x)grad.eval(x);
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

if strcmp(solver,'filterSD')
    % functions
    func = @(x)f.eval(x);
    gd = f.jacobian';
    grad = @(x)gd.eval(x);
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
