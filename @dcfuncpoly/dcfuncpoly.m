classdef dcfuncpoly < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dcfuncpoly class
    % create a dc polynomial function object
    %
    % author: yi-shuai niu
    % 2019-4
    % modified 2019-5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        f % dc function f
        g % 1se dc component
        h % 2nd dc component
        x % variables used to define a dc function
        dh=[] % dh is empty unless dh is given in constructor
        dg=[]
        df=[]
    end
    properties(Hidden)
        nvars % number of variables
    end
    properties(Dependent, SetAccess = private, GetAccess = public) % could get but not modify
        gradg % gradient of g
        gradh % gradient of h
        gradf % gradient of f
        hesseg % hessian of g
        hesseh % hessian of h
        hessef % hessian of f
    end
    methods
        % constructor
        function obj = dcfunc(x,g,h,dh)
            %%
            % obj = dcfunc
            % obj = dcfunc(x,g,h)
            % obj = dcfunc(x,g,h,dh)
            % the dc function is defined as: $f(x) = g(x) - h(x)$
            % where g and h are convex functions, and x must be a yalmip variable defined by sdpvar.
            % you can use obj.g, obj.h and obj.x to modify the class object.
            % dh is the gradient of h, if dh is not provided by user, the
            % gradient will be computed.
            % PS: Currently, for nonpolynomial function h, dh must be provided.
            if nargin == 0 % default constructor with no args
                obj.f=[];
                obj.g=[];
                obj.h=[];
                obj.x=[];
            elseif nargin == 3 % normal constructor with 3 args
                obj.f=g-h;
                obj.g=g;
                obj.h=h;
                obj.x=x;
            elseif nargin == 4 % have dh information
                obj.f=g-h;
                obj.g=g;
                obj.h=h;
                obj.x=x;
                obj.dh=dh;
            else
                error('wrong input arguments');
            end
        end
        function disp(obj)
            fprintf('----------------------------------\n')
            fprintf('DC function with %d vaiables.\n',obj.nvars);
            fprintf('----------------------------------\n')
        end
        
        % set dc objective function
        function set.f(obj,f)
            obj.f = f;
        end
        % set 1st dc component
        function set.g(obj,g)
            obj.g = g;
        end
        % set 2nd dc component
        function set.h(obj,h)
            obj.h = h;
        end
        % set variables x of the dc function
        function set.x(obj,x)
            obj.x = x;
        end
        % compute gradient of g (exact gradient)
        function gradg = get.gradg(obj)
            if isempty(obj.dg)
                gradg = jacobian(obj.g)';
            else
                gradg = obj.dg;
            end
        end
        % compute gradient of h (exact gradient)
        function gradh = get.gradh(obj)
            if isempty(obj.dh) 
                % if dh is empty, we compute gradient
                gradh = jacobian(obj.h)';
            else
                % otherwise we use user dh.
                % this is due to the fact that yalmip can not compute
                % jacobian correctly to nonlinear function.
                gradh=obj.dh;
            end
        end
        % compute gradient of f (exact gradient)
        function gradf = get.gradf(obj)
            if isempty(obj.df)
                gradf = jacobian(obj.f)';
            else
                gradf = obj.df;
            end
        end
        % compute hessian of g
        function hesseg = get.hesseg(obj)
            hesseg = hessian(obj.g,obj.x);
        end
        % compute hessian of h
        function hesseh = get.hesseh(obj)
            hesseh = hessian(obj.h,obj.x);
        end
        % compute hessian of f
        function hessef = get.hessef(obj)
            hessef = hessian(obj.f,obj.x);
        end
        % get nvars
        function nvars = get.nvars(obj)
            nvars = length(obj.x);
        end
        % evaluate g at xval
        function val = evalg(obj,xval)
            if isempty(obj.g)
                val=obj.f.eval(xval)+obj.h.eval(xval);
            else
                val=obj.g.eval(xval);
            end
        end
        % evaluate h at xval
        function val = evalh(obj,xval)
            if isempty(obj.h)
                val=obj.g.eval(xval)-obj.f.eval(xval);
            else
                val = obj.h.eval(xval);
            end
        end
        % evaluate f at xval
        function val = evalf(obj,xval)
            if isempty(obj.f)
                val=obj.g.eval(xval)-obj.h.eval(xval);
            else
                val = obj.f.eval(xval);
            end
        end
        % evaluate gradient of g at xval
        function val = evalgradg(obj,xval)
            val = obj.gradg.eval(xval);
        end
        % evaluate approximate gradient of g at xval
        function val = evalapproxgradg(obj,xval,delta)
            func=@(x)obj.evalg(x);
            if nargin < 3
                delta=0.01;
            end
            val = numgrad(func,xval,delta);
        end
        % evaluate gradient of h at xval
        function val = evalgradh(obj,xval)
            val = obj.gradh.eval(xval);
        end
        % evaluate approximate gradient of h at xval
        function val = evalapproxgradh(obj,xval,delta)
            func=@(x)obj.evalh(x);
            if nargin < 3
                delta=0.01;
            end
            val = numgrad(func,xval,delta);
        end
        % evaluate gradient of f at xval
        function val = evalgradf(obj,xval)
            val = obj.gradf.eval(xval);
        end
        % evaluate approximate gradient of f at xval
        function val = evalapproxgradf(obj,xval,delta)
            func=@(x)obj.evalf(x);
            if nargin < 3
                delta=0.01;
            end
            val = numgrad(func,xval,delta);
        end
        % evaluate hessian of g at xval
        function val = evalhesseg(obj,xval)
            val = obj.hesseg.eval(xval);
        end
        % evaluate hessian of h at xval
        function val = evalhesseh(obj,xval)
            val = obj.hesseh.eval(xval);
        end
        % evaluate hessian of f at xval
        function val = evalhessef(obj,xval)
            val = obj.hessef.eval(xval);
        end
        
        % overload operators
        % plus
        % two function must have same variables
        function newobj = plus(obj1,obj2)
            newobj = dcfuncpoly(obj1.x,obj1.g+obj2.g,obj1.h+obj2.h);
        end
        % minus
        function newobj = minus(obj1,obj2)
            newobj = dcfuncpoly(obj1.x,obj1.g+obj2.h,obj1.h+obj2.g);
        end
        % multiplication
        function newobj = mtimes(a,obj)
            if isa(a,'dcfuncpoly') %
                % under construction
            elseif isa(a,'double') % scalar times dc function
                if a>=0
                    newobj = dcfuncpoly(obj.x,a*obj.g,a*obj.h);
                else
                    newobj = dcfuncpoly(obj.x,-a*obj.h,-a*obj.g);
                end
            end
        end
        
        % dc decomposition for some special class of dc functions
        function specdcd(obj)
            if is(obj.f,'quadratic')
                % spectral dc decomposition to quadratic function
                Q=obj.hessef;
                [P,V]=eig(full(Q/2));
                v=diag(V);
                lstn=find(v<0);
                
                V1=V;
                V1(lstn,lstn)=0;
                Q1 = P*V1/P;
                g1 = obj.x'*Q1*obj.x;
                h1 = g1 - obj.f;
                obj.g = g1;
                obj.h = h1;
                fprintf('the dc function is turned into spectral dc decomposition.\n');
            else
                warning('the dc function is not quadratic!');
            end
        end
        
    end
end

