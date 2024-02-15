function constraints = problem_inequality_constraints(a, b, delta)
%PROBLEM_INEQUALITY_ CONSTRAINTS Set up inequality constraint polynomials.
%
% NB: /equality/ constraints are applied at the level of operators when
% creating the obj.Scenario object.
%

    % Check number of a matches number of b, and get Scenario
    d = numel(a);
    a.checkSameScenario(b);        
    scenario = a.Scenario;
    assert(numel(b) == d);
    
    % Make constraints: 1 - b_i^2 >= 0 for i = 1, ..., d^2
    constraints = struct;
    constraints.neg_square = 1.0 - b .* b;
    constraints.neg_square.ReadOnly = true;

    % Make constraints: δ ± i[a_j, b_k] >= 0 for j, k = 1, ..., d
    constraints.pos_comm = ...
        MTKPolynomial.InitForOverwrite(scenario, [d*d, 1]);
    constraints.neg_comm = ...
        MTKPolynomial.InitForOverwrite(scenario, [d*d, 1]);
    delta_base = MTKMonomial.InitValue(scenario, delta * ones([d, 1]));
    for idx = 1:d % Loop over a_1... a_d                
        splice = 1i * commutator(a(idx), b); % ((a(idx) .* b) - (b .* a(idx)));                
        pos_splice = delta_base + splice;
        neg_splice = delta_base - splice;
        range_idx = (((idx-1)*d)+1):(idx*d);
        constraints.pos_comm(range_idx) = pos_splice;
        constraints.neg_comm(range_idx) = neg_splice;                                
    end
    constraints.pos_comm.ReadOnly = true;
    constraints.neg_comm.ReadOnly = true;
end

