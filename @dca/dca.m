classdef dca < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dca class
    % dc programming solver
    %
    % Author: yi-shuai niu
    % 2016-3
    % 2019-3 modified
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
        convexsolver = 'ipopt' % convex subproblem solver, default cplex
        yalmipoptions; % yalmip setting
        linesearch = false % use line search for acceleration
        plotlinetype='b-s'; % if plot = 1, this option set the plot line type
        plotinnewfig=true; % if true, we flot in a new fig, otherwise, we plot in fig 1
        rho;
        boostthreshold=1e-3;
    end
    
    methods
        function obj = dca(dcp,x0)
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
        function set.yalmipoptions(obj,ops)
            % set yalmip option
            obj.yalmipoptions = ops;
        end
        function set.rho(obj,val)
            % set parameter rho for universal decomposition
            obj.rho = val;
        end
        function set.boostthreshold(obj,val)
            obj.boostthreshold = val;
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
            
            obj.iter = 0;
            xk = obj.x0;
            if strcmp(obj.convexsolver,'bpppa')==false
                 obj.yalmipoptions = sdpsettings('solver',obj.convexsolver,'verbose',0);
            end
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
                fprintf('DCA version 1.0 with yalmip \nLocal solver: %s\n',obj.convexsolver);
                if obj.linesearch
                    fprintf('* activate Armijo linesearch acceleration\n');
                end
                fprintf('------------------------------------------------------------\n');
                fprintf('Iterations | Objective values |   Delta x   |   Delta f \n');
            end
            yalmiptime=0;
            tic;
            while obj.iter < obj.maxiter
                obj.iter = obj.iter+1;
                yk = obj.dcp.F.evalgradh(xk); % compte gradient of h at xk
                if strcmp(obj.convexsolver,'bpppa')==false
                    diagnostics=solvesdp(obj.dcp.C,obj.dcp.F.g-obj.dcp.X'*yk,obj.yalmipoptions);
                else
                    xx=BPPPA(-yk/obj.rho,1);
                    assign(obj.dcp.X,xx);
                    diagnostics.problem = 0; 
                    diagnostics.yalmiptime = 0;
                end
                yalmiptime = yalmiptime + diagnostics.yalmiptime;
                % check feasibility of subproblems
                if diagnostics.problem ~=0
                    if (obj.verbose == 1)
                        fprintf('------------------------------------------------------------\n');
                    end
                    status = setstatus(toc,yalmiptime,obj.iter,2,'Problem infeasible or undounded.');
                    return;
                end
                xk1 = value(obj.dcp.X);
                fk = obj.dcp.F.evalf(xk);
                fk1 = obj.dcp.F.evalf(xk1);
                normf = abs(fk1-fk);
                normx = norm(xk1-xk);
                % accelerate with line search
                if obj.linesearch == true
                    if (normf<obj.boostthreshold && normx~=0)
                        [xacc,facc]=armijo(obj,xk1,xk1-xk);
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
                    status = setstatus(toc,yalmiptime,obj.iter,0,'Successfully solved.');
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
            status = setstatus(toc,yalmiptime,obj.iter,1,'Max interation exceed.');
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
function status = setstatus(timer,yalmiptime,iter,flag,info)
status.time = timer-yalmiptime;
status.yalt = yalmiptime;
status.iter = iter;
status.avgt = timer/iter;
status.flag = flag;
status.info = info;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% accelerator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xx = accelerator(obj)
% line search to get a better feasible solution
x0=value(obj.dcp.X);
d=-obj.dcp.F.evalgradf(x0);
a=sdpvar(1,1);
optimize(obj.dcp.C + [obj.dcp.X == x0+a*d; a>=0],obj.dcp.F.f-a)
xx=value(obj.dcp.X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% armijo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff] = armijo_old(obj,xk,xk1)
% line search to get a better feasible solution
x0=xk1;
f0=obj.dcp.F.evalf(x0);
%d=-obj.dcp.F.evalgradf(x0);
d=xk1-xk;
beta = 0.5;
s = sqrt(2)/(beta*norm(d)); % 10000
alpha=beta*s;
while (alpha>1e-8)
    xx=x0+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=f0 - ff - alpha*norm(d)^2;
    if delta >=0 && min(check(obj.dcp.C))>=-1e-7
        return;
    end
    alpha=beta*alpha;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% armijo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff] = armijo(obj,xk1,d)
% line search to get a better feasible solution
fk1=obj.dcp.F.evalf(xk1);
nd = norm(d);
beta = 0.618;
alpha = sqrt(2)/nd; % 10000
while (alpha*nd>1e-5)
    xx=xk1+alpha*d;
    ff=obj.dcp.F.evalf(xx);
    delta=fk1 - ff - alpha*nd^2;
    if delta >=-1e-8 && checkfeas(xx)
        return;
    end
    alpha=beta*alpha;
end
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