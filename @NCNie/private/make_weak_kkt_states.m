function states = make_weak_kkt_states(obj)
%MAKE_WEAK_KKT_STATES 
% Construct various states associated with objective and KKT conditions.
% We generate only the a vectors, since the problem is symmetric under
% complex conjugation (so we can 'ignore' the b vectors).
%
% Outputs the following structure, whose leaves are yalmip sdpvars:
% states
%   .sigma          [1 state]
%   .mu             [struct]
%      .neg_square  [d states]
%      .pos_comm    [d^2 states] in .a (real) and .b (imaginary) parts
%      .neg_comm    [d^2 states] in .a (real) and .b (imaginary) parts
%   .lambda         [d states]
%
% This state layout is the same for weak and essential ncKKT conditions.
%

   states = struct;

   % One state for σ.
   states.sigma = make_sigma(obj);

   % Collection of μ states for inequalities (g_i).
   states.mu = make_mu(obj);
   
   % Collection of λ states for equalities (h_i).
   states.lambda = make_lambda(obj);        
end

%% Private functions
function sigma = make_sigma(obj)
   sigma = obj.Scenario.yalmipVars;
end

function mu = make_mu(obj)
    mu = struct;
    % μ:  b_i^2 - 1 >= 0
    mu.neg_square = fill_state_cell(obj, obj.d);

    % μ: δ + i[a_j, b_k] >= 0, could be complex
    mu.pos_comm = fill_complex_state_cell(obj, obj.d * obj.d);
    
    % μ: δ - i[a_j, b_k] >= 0, could be complex 
    mu.neg_comm = fill_complex_state_cell(obj, obj.d * obj.d);
end

function lambda = make_lambda(obj)
    % λ: a_i^2 - 1 = 0
    lambda = fill_state_cell(obj, obj.d);
end

function res = fill_state_cell(obj, N)
    res = cell(1, N);
    for i = 1:N
        res{i} = obj.Scenario.yalmipVars();
    end
end

function res = fill_complex_state_cell(obj, N)
    res = cell(1, N);
    for i = 1:N
        res{i} = struct;
        [res{i}.a, res{i}.b] = obj.Scenario.yalmipVars();
    end
end