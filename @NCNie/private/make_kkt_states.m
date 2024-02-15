function states = make_kkt_states(obj)
%MAKE_KKT_STATES 
% Construct various states associated with objective and KKT conditions.
%
% Outputs the following structure, whose leaves are yalmip sdpvars:
% states
%   .sigma          [1 state]
%   .mu             [struct]
%      .max_sphere  [1 state] in .a (real) and .b (imaginary) parts
%      .min_sphere  [1 state] in .a (real) and .b (imaginary) parts
%      .comm_plus   [3 states] in .a (real) and .b (imaginary) parts
%      .comm_minus  [3 states] in .a (real) and .b (imaginary) parts
%
% If not in exterior mode, .min_sphere is not generated.
%
    
   states = struct;

   % One state for σ.
   states.sigma = make_sigma(obj);

   % Collection of μ states for various inequality constraints (g_i).
   states.mu = make_mu(obj);   
end

%% Private functions
function sigma = make_sigma(obj)
   % sigma is always real
   sigma = obj.Scenario.yalmipVars;
end

function mu = make_mu(obj)
    mu = struct;
    % μ associated with max - x_1^2 - x_2^2 - x_3^2 >= 0, could be complex
    mu.max_sphere = struct;
    [mu.max_sphere.a, mu.max_sphere.b] = obj.Scenario.yalmipVars();
    
     % μ associated with x_1^2 + x_2^2 + x_3^2 - min >= 0, could be complex
    if obj.exterior       
        mu.min_sphere = struct;
        [mu.min_sphere.a, mu.min_sphere.b] = obj.Scenario.yalmipVars();
    end
    
    % μ associated with δ + i[a_j, b_k] >= 0, could be complex
    mu.comm_plus = fill_complex_state_cell(obj, 3);
    
    % μ associated with δ - i[a_j, b_k] >= 0, could be complex 
    mu.comm_minus = fill_complex_state_cell(obj, 3);
end

function res = fill_complex_state_cell(obj, N)
    res = cell(1, N);
    for i = 1:N
        res{i} = struct;
        [res{i}.a, res{i}.b] = obj.Scenario.yalmipVars();
    end
end