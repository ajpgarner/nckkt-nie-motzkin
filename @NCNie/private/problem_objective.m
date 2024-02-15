function result = problem_objective(a, b, FC)
%PROBLEM_OBJECTIVE Set up problem's objective function polynomial.

    % Check number of a matches number of b, and get Scenario
    d = numel(a);
    a.checkSameScenario(b);            
    assert(numel(b) == d);
    scenario = a.Scenario;
    
    % Check FC dimensions
    assert(isequal(size(FC), [d+1, d+1]));   
    
    % Constant term (probably zero)
    result = MTKMonomial.InitValue(scenario, FC(1,1));
    
    % Terms in a_i or b_j
    for idx = 1:d
        if FC(1, idx+1) ~= 0
            result = result + FC(1,idx+1) * a(idx);
        end
        if FC(idx+1, 1) ~= 0
            result = result + FC(idx+1, 1) * b(idx);
        end
    end
    
    % Terms in 0.5 {a_i, b_j} (equivalently: Re(a_i b_i) )
    for idx = 1:d
        for jdx = 1:d
            if FC(idx+1, jdx+1) ~= 0
                result = result ...
                    + FC(idx+1, jdx+1) * 0.5 * anticommutator(a(idx), b(jdx));
            end
        end
    end   
    
     % Lock
    result.ReadOnly = true;
end
