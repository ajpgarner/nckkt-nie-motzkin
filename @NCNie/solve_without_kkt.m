function [result, state] = solve_without_kkt(obj, mm_level, lm_level, so)
%SOLVE_NO_KKT Solve, with no KKT operator conditions
% PARAMS:
%    mm_level - Hierarchy level of moment matrix
%    lm_level - Hierarchy level of localizing matrices
%    so_level - Maximum word length of state-optimality commutator
%               constraints (set to 0 or false to disable).
%

    % Basic validation of arguments
    assert(nargin>=3, ...
        "Arguments: moment matrix level, localizing matrix level.");
    if nargin < 4
        so = false;
    end
    if so == false
        so_level = 0;
    else
        so_level = so;
    end
    
    % MM level must be 3
    assert(mm_level >= 3, ...
           "MM Level must be at least 3 to contain all terms in the objective.");
        
    % Reset if necessary
    if obj.solve_state > 0
        obj.reset();
    end
    obj.solve_state = 1;
    
    % Report solve type
    if obj.Verbose>=1
        if so
            fprintf("\nSolving with state-optimality only, MM = %d, LM = %d, SO = %d...\n", ...
                mm_level, lm_level, so_level);
        else
            fprintf("\nSolving without ncKKT, MM = %d, LM = %d...\n", ...
                mm_level, lm_level);
        end
    end
    
    % Make moment and localizing matrices
    if obj.Verbose>=1
        fprintf("Generating matrices [MM %d, LM %d]...\n", mm_level, lm_level);    
    end
    mm_start = tic;
    [obj.mm, obj.lm] = obj.make_matrices(mm_level, lm_level);    
    if so
        gamma = obj.make_gamma_matrix(lm_level);
    end
    
    mm_time = toc(mm_start);
    if obj.Verbose>=1
        fprintf("...matrices generated in %f seconds.\n", mm_time);
    end
   
    ym_start = tic;
    if obj.Verbose>=1
        fprintf("Setting up constraints in yalmip...\n");
    end
    
    % Reset yalmip
    yalmip('clear');
    
    % Problem is symmetric under conjugation, so ignore imaginary 'b'
    sigma_a = obj.Scenario.yalmipVars();
    
    % Constrain normalization
    constraints = [sigma_a(1) == 1];
    
    % Constrain PSD moment matrix and localizing matrices
    constraints = [constraints, ...
                   obj.base_psd_conditions(sigma_a, false, obj.mm, obj.lm)];
               
    % (Optional) state-optimality constraints
    if so
        constraints = [constraints, ...
           obj.state_optimality_conditions(gamma, sigma_a, so_level)];
    end
    
    % Objective function:
    objective = obj.Objective.yalmip(sigma_a);
    
    ym_time = toc(ym_start);
    if obj.Verbose >= 1
        fprintf("...constraints generated in %f seconds.\n", ym_time);
        fprintf("Passing to solver...\n");
    end
    
    % Solver parameters
    settings = sdpsettings;
    if obj.Verbose >= 1 
        settings = sdpsettings(settings, 'verbose', true);
    else
        settings = sdpsettings(settings, 'verbose', false);
    end
    settings = sdpsettings(settings, 'solver', 'mosek');
    
    % Solve [minimization]
    ym_diagnostics = optimize(constraints, objective, settings); 
    assert(ym_diagnostics.problem ~= 1, "Problem was infeasible!");
	if ym_diagnostics.problem ~= 0
		fprintf("/!\ A problem occured while solving...\n");
		disp(ym_diagnostics);
	end

    % Get objective value
    result = value(objective);
    
    % Get solution state
    obj.soln_states= value(sigma_a);
    if nargout >= 2        
        state = obj.soln_states;
    end
        
    
    % Mark as solved
    obj.solve_state = 2;
    
    % Print solution in verbose mode
    if obj.Verbose >= 1
        fprintf("Solution mm = %d, lm = %d, so = %d: %.10g\n", ...
            mm_level, lm_level, so_level, result);
    end
end

