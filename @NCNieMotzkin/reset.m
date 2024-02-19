function reset(obj)
%RESET Resets moment objects in class, so as to provide a clean solution.
%
% Note: Any old handles to objects outside of the class will remain valid,
% but will become incompatible with the new objects created within this 
% class.
%

    % Print parameters
    if obj.Verbose >= 1
        if obj.exterior
            fprintf("Setting up scenario with min sphere = %g, max sphere = %g, delta = %g.\n",...
                    obj.min_sphere, obj.max_sphere, obj.delta);
        else
             fprintf("Setting up scenario with max sphere = %g, delta = %g.\n",...
                     obj.max_sphere, obj.delta);
        end
    end
    
    % Set up scenario with no (operator level) equality constraints:
    obj.Scenario = AlgebraicScenario(["x1", "x2", "x3"]);    
    obj.Scenario.System;

    % Bind operator references
    [obj.x1, obj.x2, obj.x3] = obj.Scenario.getAll();
        
    % Useful powers
    obj.x1_pow2 = obj.x1^2;
    obj.x1_pow4 = obj.x1^4;    
    obj.x2_pow2 = obj.x2^2;
    obj.x2_pow4 = obj.x2^4;
    obj.x3_pow2 = obj.x3^2;
    obj.x3_pow4 = obj.x3^4;
    obj.x3_pow6 = obj.x3^6;
    
    icomm_x1x2 = 1.0i * (obj.x1 * obj.x2 - obj.x2 * obj.x1);
    icomm_x2x3 = 1.0i * (obj.x2 * obj.x3 - obj.x3 * obj.x2);
    icomm_x3x1 = 1.0i * (obj.x3 * obj.x1 - obj.x1 * obj.x3);
    
    % Make symmetrized objective polynomial    
    obj.Objective = 0.5*anticommutator(obj.x1_pow4, obj.x2_pow2) ...
                  + 0.5*anticommutator(obj.x1_pow2, obj.x2_pow4) ...
                  + obj.x3_pow6 ...
                  - 1.5 * obj.x1_pow2 * obj.x2_pow2 * obj.x3_pow2 ...
                  - 1.5 * obj.x3_pow2 * obj.x2_pow2 * obj.x1_pow2;
              
    % Make inequality constraints
    obj.Constraints = struct;
    obj.Constraints.max_sphere = obj.max_sphere - obj.x1_pow2 - obj.x2_pow2 - obj.x3_pow2;
    
    if (obj.exterior)
        obj.Constraints.min_sphere = obj.x1_pow2 + obj.x2_pow2 + obj.x3_pow2 - obj.min_sphere;
    end
        
    obj.Constraints.comm_plus = cell(3, 1);
    obj.Constraints.comm_plus{1} = obj.delta + icomm_x2x3;
    obj.Constraints.comm_plus{2} = obj.delta + icomm_x3x1;
    obj.Constraints.comm_plus{3} = obj.delta + icomm_x1x2;
    
    obj.Constraints.comm_minus = cell(3, 1);
    obj.Constraints.comm_minus{1} = obj.delta - icomm_x2x3;
    obj.Constraints.comm_minus{2} = obj.delta - icomm_x3x1;
    obj.Constraints.comm_minus{3} = obj.delta - icomm_x1x2;
            
    % Mark as unsolved and clear solution information
    obj.solve_state = 0;
    obj.soln_states = false;
    obj.mm = false;
    obj.lm = false;
    obj.gamma = false;

end