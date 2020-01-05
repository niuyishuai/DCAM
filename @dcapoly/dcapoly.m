classdef dcapoly < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dca class for polynomial optimization (without Yalmip)
    % dc programming solver
    %
    % Author: yi-shuai niu
    % 2019-4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        dcp % dc program
        x0 % initial point
    end
    properties(GetAccess = public,SetAccess = private)         % read only
        fopt = inf % objective value
        xopt = [] % optimal solution
        iter = 0 % iterations of dca
    end
    properties
        tolf = 1e-6 % tolerence for objective function
        tolx = 1e-6 % tolerence for iterative point
        maxiter = 10000 % max iterations for dca
        plot = 0  %1: draw iterations of dca, 0: otherwise
        verbose = 1  %1: display iterations of dca, 0: otherwise
        convexsolver = 'knitro' % convex subproblem solver, default cplex
        linesearch = false % use line search for acceleration
        plotlinetype='b-s'; % if plot = 1, this option set the plot line type
        plotinnewfig=true; % if true, we flot in a new fig, otherwise, we plot in fig 1
        rho;
        boostthreshold=1e-2;
        approxgrad=false; % using approximate gradient (for large-scale cases)
        stepsize=inf; % initial step size for armijo rule
        usearmijo_adaptive=true; % use armijo adaptive rule
        armojo_adaptive_method=0; % 1: increase only once, 0: increase many times
    end
    
    methods
        function obj = dcapoly(dcp,x0)
            % dca constructor
            % obj = dca(dcp,x0)
            % where dcp is a dc program object
            % x0 is an initial point. It will be a random point if x0 is not given.
            if nargin==1
                obj.dcp = dcp;
                x0 = rand(size(dcp.X));
                obj.x0 = x0(:);
            elseif nargin==2
                obj.dcp = dcp;
                obj.x0 = x0(:);
            else
                error('wrong input arguments.');
            end
        end
        function xopt = get.xopt(obj)
            % get optimal solution
            xopt = obj.xopt;
        end
        function fopt = get.fopt(obj)
            % get objective value
            fopt = obj.fopt;
        end
        function set.tolf(obj,val)
            % set tolerence of objective function
            obj.tolf = val;
        end
        function set.tolx(obj,val)
            % set tolerence of iterative point
            obj.tolx = val;
        end
        function set.maxiter(obj,val)
            % set max iterations of dca
            obj.maxiter = val;
        end
        function set.plot(obj,val)
            % set plotting option of dca
            obj.plot = val;
        end
        function set.verbose(obj,val)
            % set verbose option of dca
            obj.verbose = val;
        end
        function set.convexsolver(obj,val)
            % set convex subproblem solver (used for yalmip)
            obj.convexsolver = val;
        end
        function set.linesearch(obj,yn)
            % set line search
            obj.linesearch = yn;
        end
        function set.rho(obj,val)
            % set parameter rho for universal decomposition
            obj.rho = val;
        end
        function set.boostthreshold(obj,val)
            obj.boostthreshold = val;
        end
        function set.approxgrad(obj,val)
            obj.approxgrad = val;
        end
        function set.usearmijo_adaptive(obj,val)
            obj.usearmijo_adaptive = val;
        end
        function set.armojo_adaptive_method(obj,val)
            obj.armojo_adaptive_method = val;
        end
        % dca algorithm for solving dcp with starting point x0
        function status = optimize(obj)
            % dca optimizer
            % status.flag : 0 dca converges with tolx or tolf
            %               1 maxiter exceed
            %               2 problem infeasible or unbounded
            % status.info : solution informations.
            % status.iter : number of iterations of dca.
            % status.time : cpu time for dca (sec.)
            % status.yalt : yalmip time (sec.)
            % status.avgt : average time for each iteration (sec.)
            
            n=numel(obj.x0);
            options = optimset('Algorithm', 'interior-point', 'Display','off');
            
            obj.iter = 0;
            xk = obj.x0;
            % plotting if actived
            if (obj.plot==1)
                if (obj.plotinnewfig)
                    figure
                else
                    figure(1);
                end
            end
            % display of actived
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
                fprintf('DCA version 1.0 beta \nLocal solver: %s\n',obj.convexsolver);
                if obj.linesearch
                    fprintf('* activate Armijo linesearch acceleration\n');
                end
                fprintf('------------------------------------------------------------\n');
                fprintf('Iterations | Objective values |   Delta x   |   Delta f \n');
            end
            cputime=tic;
            
            consequtive_count=0;
            while obj.iter < obj.maxiter
                obj.iter = obj.iter+1;
                if obj.approxgrad
                    yk = obj.dcp.F.evalapproxgradh(xk); % compte numerical gradient of h at xk
                else
                    yk = obj.dcp.F.evalgradh(xk); % compte gradient of h at xk
                end
                
                %*********************************
                if strcmpi(obj.convexsolver,'bpppa')
                    diagnostics.x=BPPPA(-yk/obj.rho,1);
                    diagnostics.problem = 0;
                else
                    diagnostics=callconvexsolver_bis(obj.dcp.F.g,-yk,obj.dcp.C, xk, obj.convexsolver,options);
                end
                % check feasibility of subproblems
                if diagnostics.problem ~=0
                    if (obj.verbose == 1)
                        fprintf('------------------------------------------------------------\n');
                    end
                    status = setstatus(toc(cputime),obj.iter,2,'Problem infeasible or undounded.');
                    return;
                end
                xk1 = diagnostics.x;
                fk = obj.dcp.F.evalf(xk);
                fk1 = obj.dcp.F.evalf(xk1);
                normf = abs(fk1-fk);
                normx = norm(xk1-xk);
                % accelerate with line search
                if obj.linesearch == true
                    if (normf<obj.boostthreshold && normx~=0)
                        if obj.usearmijo_adaptive
                            if obj.armojo_adaptive_method==0 % forwardtracking many times
                                 [xacc,facc,obj.stepsize]=armijo_adaptive(obj,xk1,xk1-xk,obj.stepsize); % increase many times
                            elseif obj.armojo_adaptive_method==1 % forwardtracking once
                                [xacc,facc,obj.stepsize]=armijo_adaptive_once(obj,xk1,xk1-xk,obj.stepsize); % increase only once
                            else
                                [xacc,facc,obj.stepsize,consequtive_count]=armijo_adaptive_s(obj,xk1,xk1-xk,obj.stepsize,consequtive_count); % increase only if two cons
                            end
                        else
                            % armijo
                            [xacc,facc]=armijo(obj,xk1,xk1-xk);
                        end
                        %if min(check(obj.dcp.C))>-1e-7
                        if obj.verbose == 1
                            fprintf('accelerated: reduced %17.3e  moved %17.3e \n',facc-fk1,norm(xacc-xk1));
                        end
                        fk1 = facc;
                        xk1 = xacc;
                    end
                end
                % compute errors
                normx = norm(xk1-xk);
                normf = abs(fk1-fk);
                % display iterations of dca if actived
                if (obj.verbose == 1)
                    fprintf('%5d %19.5e %17.3e %13.3e\n',obj.iter,fk1,normx,normf);
                end
                % plotting if actived
                if (obj.plot==1)
                    myplotf(fk,fk1,obj.iter,obj.plotlinetype);
                end
                % check stopping
                if (normx < obj.tolx || normf < obj.tolf)
                    if (obj.verbose == 1)
                        fprintf('------------------------------------------------------------\n');
                    end
                    obj.fopt = fk1;
                    obj.xopt = xk1;
                    status = setstatus(toc(cputime),obj.iter,0,'Successfully solved.');
                    return;
                end
                
                xk = xk1;
            end
            % maxiter exceed
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
            end
            obj.fopt = fk1;
            obj.xopt = xk1;
            status = setstatus(toc(cputime),obj.iter,1,'Max interation exceed.');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myplotf(fk,fk1,iter,plotline)
hold on;
if iter == 1
    title('DCA iterations');
    xlabel('Iterations');
    ylabel('Objectives');
end
if iter>1
    plot([iter-1,iter], [fk,fk1],plotline);
    drawnow
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting dca solution status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = setstatus(timer,iter,flag,info)
status.time = timer;
status.iter = iter;
status.avgt = timer/iter;
status.flag = flag;
status.info = info;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% armijo
% classical armijo rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff] = armijo(obj,xk1,d)
fk1=obj.dcp.F.evalf(xk1);
nd = norm(d);
sigma=1e-3;
beta = 0.3;
alpha = sqrt(2)/nd; % 10000
while (alpha*nd>1e-5)
    xx=xk1+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=fk1 - ff - sigma*alpha^2*nd^2;
    if delta >0 && checkfeas(xx)
        return;
    end
    alpha=beta*alpha;
end
xx=xk1;
ff=fk1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive armijo (method=0)
% increase once ()many times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff,cur_alpha] = armijo_adaptive_once(obj,xk1,d,cur_alpha)
% line search to get a better feasible solution
nd = norm(d);
alpha = min(cur_alpha,sqrt(2)/nd);
fk1=obj.dcp.F.evalf(xk1);
beta = 0.3;
backtrack=false;
forwardtrack=false;
while (alpha*nd>1e-5)
    xx=xk1+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=fk1 - ff; % - alpha^2*nd^2;
    if delta < 1e-8 || ~checkfeas(xx)
        alpha=beta*alpha; % back-tracking
        backtrack=true;
    else
        if backtrack==true || forwardtrack==true
            cur_alpha=alpha;
            return;
        else
            alpha=alpha/beta; % forward-tracking
            forwardtrack=true;
        end
    end
end
xx=xk1;
ff=fk1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive armijo (method=0)
% increase many times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff,cur_alpha] = armijo_adaptive(obj,xk1,d,cur_alpha)
% line search to get a better feasible solution
nd = norm(d);
alpha = min(cur_alpha,sqrt(2)/nd);
fk1=obj.dcp.F.evalf(xk1);
beta = 0.3;
backtrack=false;
while (alpha*nd>1e-5)
    xx=xk1+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=fk1 - ff; % - alpha^2*nd^2;
    if delta < 1e-8 || ~checkfeas(xx)
        alpha=beta*alpha; % back-tracking
        backtrack=true;
    else
        if backtrack==true
            cur_alpha=alpha;
            return;
        else
            alpha=alpha/beta; % forward-tracking
        end
    end
end
xx=xk1;
ff=fk1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive armijo (method=0)
% increase only if two consequtive iterate with same step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff,cur_alpha,consequtive_count] = armijo_adaptive_s(obj,xk1,d,cur_alpha,consequtive_count)
% line search to get a better feasible solution
nd = norm(d);
alpha = min(cur_alpha,sqrt(2)/nd);
fk1=obj.dcp.F.evalf(xk1);
beta = 0.3;
backtrack=false;
forwardtrack=false;
while (alpha*nd>1e-5)
    xx=xk1+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=fk1 - ff; % - alpha^2*nd^2;
    if delta < 1e-8 || ~checkfeas(xx)
        alpha=beta*alpha; % back-tracking
        backtrack=true;
    else
        if backtrack==true || forwardtrack==true
            cur_alpha=alpha;
            consequtive_count=0; % reset cons
            return;
        else
            if consequtive_count==2
                alpha=alpha/beta; % forward-tracking
                forwardtrack=true;
                consequtive_count=0; % reset consequtive_count
            else
                consequtive_count=consequtive_count+1;
                cur_alpha=alpha;
                return;
            end
        end
    end
end
xx=xk1;
ff=fk1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive armijo (method=1)
% increase only once
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff,alpha] = armijo_adaptive_v1(obj,xk1,d,alpha)
% line search to get a better feasible solution
fk1=obj.dcp.F.evalf(xk1);
nd = norm(d);
beta = 0.5;
if alpha==inf
    % 初始化alpha
    alpha = sqrt(2)/nd; % 10000 
end
trialtimes=0;
increasetimes=0;
while (alpha*nd>1e-5)
    xx=xk1+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=fk1 - ff - alpha^2*nd^2;
    if delta >=-1e-8 && checkfeas(xx) % check if better feasible solution
        if increasetimes == 0 && trialtimes > 0 % if never increase step size and we have tried more than once
            return;
        elseif increasetimes == 0 % if never increased
            alpha=alpha/beta; % increase alpha
            increasetimes=increasetimes+1;
        else
            return;
        end
    elseif increasetimes >0
        alpha=beta*alpha;
        xx=xk1+alpha*d;
        ff=obj.dcp.F.evalf(xx);
        return;
    else
        alpha=beta*alpha;
    end
    trialtimes=trialtimes+1;
end
alpha = inf; % reset alpha
xx=xk1;
ff=fk1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive armijo
% increase many times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff,alpha] = armijo_adaptive_v2(obj,xk1,d,alpha)
% line search to get a better feasible solution
fk1=obj.dcp.F.evalf(xk1);
nd = norm(d);
beta = 0.1;
%alpha = min(alpha,sqrt(2)/nd);
if alpha==inf
    % 初始化alpha
    alpha = sqrt(2)/nd; % 10000 
end
trialtimes=0;
increasetimes=0;
while (alpha*nd>1e-5)
    xx=xk1+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=fk1 - ff - alpha*nd^2;
    if delta >=-1e-8 && checkfeas(xx) % check if better feasible solution
        if increasetimes == 0 && trialtimes > 0 % if never increase step size and we have tried more than once
            return;
        else % if alread increased or never increase but feasible at inital size
            alpha=alpha/beta; % increase alpha
            increasetimes=increasetimes+1;
        end
    elseif increasetimes >0
        alpha=beta*alpha;
        xx=xk1+alpha*d;
        ff=obj.dcp.F.evalf(xx);
        return;
    else
        alpha=beta*alpha;
    end
    trialtimes=trialtimes+1;
end
alpha = inf; % reset alpha
xx=xk1;
ff=fk1;
end

function r=checkfeas(x)
%min(check(obj.dcp.C))>=-1e-7
if abs(sum(x)-1) <=1e-7 && min(x) >= -1e-7
    r=true;
else
    r=false;
end
end