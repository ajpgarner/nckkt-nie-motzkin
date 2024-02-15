function reset(obj)
%RESET Resets moment objects in class, so as to provide a clean solution.
%
% Note: Any old handles to objects outside of the class will remain valid,
% but will become incompatible with the new objects created within this 
% class.
%

    % Set up scenario with no (operator level) equality constraints:
    obj.Scenario = AlgebraicScenario(["x1", "x2", "x3"]);    
    obj.Scenario.System;

    % Bind operator references
    obj.x = obj.Scenario.getAll();
    [x1, x2, x3] = obj.Scenario.getAll();
    
    % Useful powers
    x1_pow2 = x1*x1;
    x1_pow4 = x1*x1*x1*x1;
    x2_pow2 = x2*x2;
    x2_pow4 = x2*x2*x2*x2;
    x3_pow2 = x3*x3;
    x3_pow6 = x3*x3*x3*x3*x3*x3;
    
    % Make objective polynomial    
    obj.Objective = x1_pow4 * x2_pow2 ...
                  + x1_pow2 * x2_pow4 ...
                  + x3_pow6 ...
                  - 3.0 * x1_pow2 * x2_pow2 * x3_pow2;
    
    
    % Make inequality constraints
    obj.Constraints = struct;
    obj.Constraints.sphere = 1 - x1_pow2 - x2_pow2 - x3_pow2;
    obj.Constraints.comm_plus = cell(3, 1);
    obj.Constraints.comm_plus{1} = obj.delta + 1i*(x1*x2 - x2*x1);
    obj.Constraints.comm_plus{2} = obj.delta + 1i*(x3*x1 - x1*x3);
    obj.Constraints.comm_plus{3} = obj.delta + 1i*(x2*x3 - x3*x2);
    obj.Constraints.comm_minus = cell(3, 1);
    obj.Constraints.comm_minus{1} = obj.delta - 1i*(x1*x2 - x2*x1);
    obj.Constraints.comm_minus{2} = obj.delta - 1i*(x3*x1 - x1*x3);
    obj.Constraints.comm_minus{3} = obj.delta - 1i*(x2*x3 - x3*x2);
        
    % Mark as unsolved
    obj.solve_state = 0;
end

