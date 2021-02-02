% A test example for solving a quadratic optimization problem using DCAM
% min f(x) = g(x) - h(x)
% s.t. sum(x) == 2; -10 <= x <= 5.
%
% The quadratic objective function f(x) = x'*Q*x 
% where Q is randomly generated and may not be positive semi-definit.
%
% NB: Yalmip is required, options of yalmip supported convex quadratic solvers: 
% gurobi|cplex|quadprog|ipopt|knitro|filterSD
%% problem data
n=5; % number of variables
Q1=rand(n,n); 
Q1=Q1'*Q1;
Q2=rand(n,n);
Q2=Q2'*Q2;
Q = Q1 - Q2;

%% define a DC programming model
% define a dc function through Yalmip to choose local solver
x=sdpvar(n,1); % define variables through yalmip
% dc decomposition of a quadratic function as g-h
g=x'*Q1*x;
h=x'*Q2*x;

% define a dc function
F=dcfunc(x,g,h)

% define a convex constraint
C = [-10<=x<=5; sum(x)==2];

% create a dc problem
mydcp_sp=dcp(F,C)

%% solve problem using classical DCA
x0=ones(n,1);
% you can set initial point x0 in second argument of dca, or use default random initialization
mydca = dca(mydcp_sp,x0);
% your can set parameters for dca or using default parameters
mydca.plot=1; % 1 for plot, 0 non
mydca.plotlinetype='b-o'; % set line style in plot
mydca.verbose = 1; % 1 to show verboase, 0 silence
mydca.tolf = 1e-6; % you may use default parameter
mydca.tolx = sqrt(mydca.tolf); 
% in general, you need to choose a suitable convex optimization solver installed on your pc
% for convex quadratic optimization, you can choose:
% gurobi|cplex|quadprog|ipopt|knitro|filterSD
% if you have not installed any external solver, you can use 'quadprog' of MATLAB
mydca.convexsolver='quadprog'; 
% optimization using dca
status = mydca.optimize();

% show result
mydca.xopt
mydca.fopt
mydca.iter
status

%% solve problem using boosted-DCA
mydca.linesearch=true;
mydca.plotinnewfig=false; % draw in the same figure if false
mydca.plotlinetype='r-s';
status = mydca.optimize();

%% reformulate problem using spectral DC decomposition and solve by classical DCA 
% spectral dc decomposition of F
F.specdcd;
% create dc problem
mydcp_sp=dcp(F,C);

% optimize using dca
mydca_sp = dca(mydcp_sp,x0)
mydca_sp.plot=1; % 1 for plot, 0 non
mydca_sp.plotinnewfig=false;
mydca_sp.plotlinetype='m-d';
mydca_sp.verbose = 1; % 1 to show verboase, 0 silence
mydca.tolf = 1e-6; % you may use default parameter
mydca_sp.optimize()

% show result
mydca_sp.xopt;
mydca_sp.fopt