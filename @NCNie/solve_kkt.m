function [result, state_values] = ...
    solve_kkt(obj, mm_level, lm_level, kkt_level, epsilon, so_level)
%SOLVE_KKT Solve, subject to weak and essential KKT conditions
% PARAMS:
%    mm_level - Hierarchy level of moment matrix
%    lm_level - Hierarchy level of localizing matrices
%   kkt_level - Maximum degree of monomial variates
%     epsilon - Positive real factor in mu(MM + e\id)>=0 constraints.
%    so_level - Maximum monomial length in state-optimality conditions.
%

    % Validate parameters
    assert(nargin>=5, ...
        "Arguments: moment matrix level, localizing matrix level, KKT level, epsilon.");
    if nargin < 6
        so_level = 0;
    end
    
    % Reduce to solve_without_kkt special case if kkt_level == 0:
    if kkt_level <= 0
        [result, state_values] = ...
            obj.solve_without_kkt(mm_level, lm_level, so_level);
        return;
    end
    
    % Reset if necessary
    if obj.solve_state > 0
        obj.reset();
    end
    obj.solve_state = 1;
    
    % Report solve type
    if obj.Verbose>=1
        if so_level > 0
            fprintf("\nSolving with ncKKT and state-optimality, MM = %d, LM = %d, KKT = %d, epsilon = %f, SO = %d...\n", ...
                mm_level, lm_level, kkt_level, epsilon, so_level);
        else
            fprintf("\nSolving with ncKKT, MM = %d, LM = %d, KKT = %d, epsilon = %f...\n", ...
                mm_level, lm_level, kkt_level, epsilon);
        end
    end
    
    % Make moment and localizing matrices
    if obj.Verbose>=1
        fprintf("Generating matrices...\n");
    end
    
    mm_start = tic;
    [mm, lm] = obj.make_matrices(mm_level, lm_level);    
    if so_level
        gamma = obj.make_gamma_matrix(lm_level);
    end
    mm_time = toc(mm_start);
    if obj.Verbose>=1
        fprintf("...matrices generated in %f seconds.\n", mm_time);
    end
   
    % Reset yalmip
    yalmip('clear');
    
    % Create states
    state_start = tic;
    if obj.Verbose>=1
        fprintf("Initiating states...\n");
    end    
    states = make_kkt_states(obj);
    state_time = toc(state_start);
    if obj.Verbose>=1
        fprintf("...states initiated in %f seconds.\n", state_time);
    end
    
    constraint_start = tic;
    if obj.Verbose>=1
        fprintf("Generating constraints...\n");
    end    
    
    % Normalization constraint on sigma   
    constraints = [states.sigma(1) == 1];
    
    % Positivity constraints on MM and LM (all states)
    constraints = [constraints, ...
        kkt_psd_conditions(obj, states, mm, lm, epsilon)];
            
    % mu_i(g_i) and optimality constraints
    constraints = [constraints, ...
                   kkt_state_conditions(obj, states, kkt_level)];
                             
    % (Optional) state-optimality constraints
    if so_level
        constraints = [constraints, ...
           state_optimality_conditions(obj, gamma, states.sigma, kkt_level)];
    end
               
    % In verbose mode, print a summary of constraints
    constraint_time = toc(constraint_start);
    if obj.Verbose >= 1
        cs = size(constraints);
        fprintf("...generated %d constraints in %f seconds.\n", ...
                cs(1), constraint_time);
    end
    
    % Objective function
    objective = obj.Objective.yalmip(states.sigma);
    
    % Solver parameters
    settings = sdpsettings;
    if obj.Verbose >= 1 
        settings = sdpsettings(settings, 'verbose', true);
    else
        settings = sdpsettings(settings, 'verbose', false);
    end
    settings = sdpsettings(settings, 'solver', 'mosek');
    
    % Solve
    if obj.Verbose >= 1
        fprintf("Passing to solver...\n");
    end
    
    ym_diagnostics = optimize(constraints, objective, settings);
    assert(ym_diagnostics.problem ~= 1, "Problem was infeasible!");
    if ym_diagnostics.problem ~= 0
        fprintf("/!\ A problem occured while solving...\n");
        disp(ym_diagnostics);
    end
    
    % Get objective value
    result = value(objective);

    % Get state values (nb: same layout as weak)
    obj.soln_states = evaluate_kkt_states(obj, states); 
    if nargout >=2 
        state_values = obj.soln_states;
    end
    
    % Mark as solved
    obj.solve_state = 2;

    % Print solution in verbose mode
    if obj.Verbose >= 1
        fprintf("Solution mm = %d, lm = %d, kkt = %d, e= %d, so = %d: %.10g\n", ...
            mm_level, lm_level, kkt_level, epsilon, so_level, result);
    end
end

