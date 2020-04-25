classdef dcppoly < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dcp class
    % create a dc programming problem object
    % min {F(X): X\in C}
    %
    % author: yi-shuai niu
    % 2017-3-22
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        F % objectif function
        C % constraints
        X % variables
    end
    
    properties(Dependent,Access=private)
        nvars % number of variables
        ncons % number of constraints
    end
    
    methods
        function obj = dcppoly(F,C)
            % dc program constructor
            % obj = dcp(F,C)
            % where F is the dc objective function (dcfunc class)
            % C is convex constraint (yalmip format)
            if nargin == 0
                obj.F = dcfuncpoly;
                obj.C = [];
                obj.X = [];
            else
                obj.F = F;
                obj.C = C;
                obj.X = F.x;
            end
        end
        
        function disp(obj)
            % display function of dcp
            fprintf('----------------------------------\n')
            fprintf('DC programming problem with %d variable(s) and %d constraint(s).\n',obj.nvars, obj.ncons);
            fprintf('----------------------------------\n')
        end
        
        function nvars = get.nvars(obj)
            % get number of variables
            nvars = obj.F.nvars;
        end
        function ncons = get.ncons(obj)
            % get number of constraints
            ncons = length(obj.C);
        end
    end
end