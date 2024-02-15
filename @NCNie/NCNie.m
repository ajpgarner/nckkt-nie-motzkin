classdef NCNie < handle
%NCNie Example based on example 5.3 and 5.4 of Nie [arXiv:1006.2418]
%
% To make and solve an example, first invoke the constructor of this class
% to set up the problem, and then the appropriate solve_ function to
% generate and solve an SDP relaxation of its solution.
%
    
    %% Scenario parameters
    properties(SetAccess=private, GetAccess=public)        
        delta           % Parameter in commutation constraint.
        exterior
        min_sphere      % Inner sphere bound (implictly 0 if not exterior)
        max_sphere      % Outer sphere bound.
    end
    
    %% Moment objects
    properties(SetAccess=private, GetAccess=public)
        Scenario        % Algebraic scenario (Moment handle).        
        x1              % Operator x1
        x2              % Operator x2
        x3              % Operator x3
        x1_pow2         % Utility monomial: x1^2
        x1_pow4         % Utility monomial: x1^4
    	x2_pow2         % Utility monomial: x2^2
    	x2_pow4         % Utility monomial: x2^4
    	x3_pow2         % Utility monomial: x3^2
        x3_pow4         % Utility monomial: x3^4
    	x3_pow6         % Utility monomial: x3^6
        Objective       % Objective function polynomial.
        Constraints     % Inequality constraint polynomials.
    end
    
    %% Read/write parameters
    properties(Access=public)
        Verbose = 0 % Verbosity level (0, 1 or 2).
    end
        
    %% Solution parameters
    properties(SetAccess=private, GetAccess=public)
        mm; % Moment matrix
        lm; % Localizing matrices.
    end
    
    %% Private parameters
    properties(Access=private)
        solve_state = 0; % 0 = no solve, 1 = started; 2 = completed.
    end
       
    
    %% Constructor
    methods        
        function obj = NCNie(delta, min_sphere, max_sphere, verbose)
        %MFCQEXAMPLE Construct an example problem.
        %
        % PARAMS:
        %            delta - The limit to non-commutation.
        %
                     
            % Parse delta parameter            
            assert(nargin>=1 && isscalar(delta) && delta >= 0)
            obj.delta = double(delta);
            
            if (nargin < 3)
                max_sphere = 1;
            end
            
            % Validate minimum sphere (0 for unconstrained)
            if (nargin < 2)
                min_sphere = 0;                
            end
            assert(isscalar(min_sphere) && min_sphere >= 0, ...
                "Interior sphere must have non-negative radius");
            obj.exterior = (min_sphere ~= 0);
            obj.min_sphere = min_sphere;

            % Validate maximum sphere
            assert(isscalar(max_sphere) && (max_sphere > min_sphere), ...
                        "Max sphere must be greater than min sphere.");            
            obj.max_sphere = max_sphere;
            
            % Parse verbosity setting
            if nargin >= 4
                obj.Verbose = verbose;
            else
                obj.Verbose = 0;
            end
            
            % Set-up problem
            obj.reset();
        end
    end  
end

