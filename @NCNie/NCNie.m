classdef NCNie < handle
%MFCQEXAMPLE Class for solving the simple example in arxiv:2311.18707.
%
% This problem is designed by construction to meet the non-commutative 
% Mangasarian-Fromovitz Constraint Qualification (MFCQ).
%
% To make and solve an example, first invoke the constructor of this class
% to set up the problem, and then the appropriate solve_ function to
% generate and solve an SDP relaxation of its solution.
%
    
    %% Read-only parameters
    properties(SetAccess=private, GetAccess=public)        
        delta           % Parameter in commutation constraint.
        Scenario        % Algebraic scenario (Moment handle).
        x               % Operators x1, x2 and x3.
        Constraints     % Inequality constraint polynomials.
        Objective       % Objective function polynomial.
    end
    
    %% Read/write parameters
    properties(Access=public)
        Verbose = 0 % Verbosity level (0, 1 or 2).
    end
    
    %% Private parameters
    properties(Access=private)
        solve_state = 0; % 0 = no solve, 1 = started; 2 = completed.
    end
       
    
    %% Constructor
    methods        
        function obj = NCNie(delta, verbose)
        %MFCQEXAMPLE Construct an example problem.
        %
        % PARAMS:
        %            delta - The limit to non-commutation.
        %
                     
            % Parse delta parameter            
            assert(nargin>=1 && isscalar(delta)  && isscalar(delta))
            obj.delta = double(delta);            
         
            % Parse verbosity setting
            if nargin >= 2
                obj.Verbose = verbose;
            else
                obj.Verbose = 0;
            end
            
            % Set-up problem
            obj.reset();
        end
    end  
end

