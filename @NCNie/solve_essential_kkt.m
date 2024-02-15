function [result, state_values] = ...
    solve_essential_kkt(obj, mm_level, lm_level, kkt_level, epsilon, so)
%SOLVE_ESSENTIAL_KKT Solve, subject to weak KKT conditions
% PARAMS:
%    mm_level - Hierarchy level of moment matrix
%    lm_level - Hierarchy level of localizing matrices
%   kkt_level - Maximum degree of monomial variates
%     epsilon - Positive real factor in mu(MM + e\id)>=0 constraints.
%          so - True to enable state-optimality conditions.
%

    % Validate parameters
    assert(nargin>=5, ...
        "Arguments: moment matrix level, localizing matrix level, KKT level, epsilon.");
    if nargin < 6
        so = false;
    end
    
    % Reset if necessary
    if obj.solve_state > 0
        obj.reset();
    end
    obj.solve_state = 1;
    
    % Report solve type
    if obj.Verbose>=1
        if so
            fprintf("\nSolving essential ncKKT with state-optimality, MM = %d, LM = %d, KKT = %d, epsilon = %f...\n", ...
                mm_level, lm_level, kkt_level, epsilon);
        else
            fprintf("\nSolving essential ncKKT, MM = %d, LM = %d, KKT = %d, epsilon = %f...\n", ...
                mm_level, lm_level, kkt_level, epsilon);
        end
    end
    
    % Make moment and localizing matrices
    if obj.Verbose>=1
        fprintf("Generating matrices...\n");
    end
    
    mm_start = tic;
    [mm, lm] = obj.make_matrices(mm_level, lm_level);
    [mm_eps, lm_eps] = obj.make_epsilon_matrices(mm, lm, epsilon);    
    if so
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
    states = make_weak_kkt_states(obj); % Essential same as weak!
    state_time = toc(state_start);
    if obj.Verbose>=1
        fprintf("...states initiated in %f seconds.\n", state_time);
    end
    
    constraint_start = tic;
    if obj.Verbose>=1
        fprintf("Generating constraints...\n");
    end    
    
    % Normalization constraint    
    constraints = [states.sigma(1) == 1];
    
    % Positivity constraints on MM and LM (all states)
    constraints = [constraints, ...
                   essential_kkt_psd_conditions(obj, states, mm, lm, ...
                                                mm_eps, lm_eps)];
            
    % Normed-KKT equality constraints
    constraints = [constraints, ...
                   essential_kkt_state_conditions(obj, states, kkt_level)];
                             
    % (Optional) state-optimality constraints
    if so
        constraints = [constraints, ...
           obj.state_optimality_conditions(gamma, states.sigma, kkt_level)];
    end
               
    % In verbose mode, print a summary of constraints
    constraint_time = toc(constraint_start);
    if obj.Verbose >= 1
        cs = size(constraints);
        fprintf("...generated %d constraints in %f seconds.\n", ...
                cs(1), constraint_time);
        if obj.Verbose >= 2 
            % If very verbose, also display yalmip constraint table:
            fprintf("\n");
            constraints %#ok<NOPRT>
        end
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
    ym_diagnostics = optimize(constraints, objective, settings); 
    assert(ym_diagnostics.problem ~= 1, "Problem was infeasible!");
         
    % Get objective value
    result = value(objective);

    % Get state values (nb: same layout as weak)
    state_values = evaluate_weak_kkt_states(obj, states); 
    
    % Mark as solved
    obj.solve_state = 2;

    % Print solution in verbose mode
    if obj.Verbose >= 1
        fprintf("Solution: %.10g\n", result);
    end
end

