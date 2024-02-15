function val_states = evaluate_normed_kkt_states(obj, ym_states)
%EVALUATE_NORMED_KKT_STATES Get actual numbers associated with solution.

    val_states = struct;
    val_states.sigma = value(ym_states.sigma);
    
    % Collection of mu states
    val_states.mu  = evaluate_mu(obj, ym_states.mu);
    
    % Collection of lambda+ states
    val_states.lambda_plus  = evaluate_lambda(obj, ym_states.lambda_plus);
        
    % Collection of lambda- states
    val_states.lambda_minus = evaluate_lambda(obj, ym_states.lambda_minus);
    
end

%% Private functions
function mu = evaluate_mu(obj, ym_states)
    mu = struct;
    mu.neg_square = evaluate_states(obj.d, ym_states.neg_square);
    mu.pos_comm = evaluate_complex_states(obj.d * obj.d, ym_states.pos_comm);
    mu.neg_comm = evaluate_complex_states(obj.d * obj.d, ym_states.neg_comm);
end

function lambda = evaluate_lambda(obj, ym_states)
    lambda = evaluate_states(obj.d, ym_states);    
end

function res = evaluate_states(N, ym_states)
    res = cell(1, N);
    for i = 1:N
        res{i} = value(ym_states{i});
    end
end

function res = evaluate_complex_states(N, ym_states)
    res = cell(1, N);
    for i = 1:N
        res{i} = struct;
        res{i}.a = value(ym_states{i}.a);
        res{i}.b = value(ym_states{i}.b);
    end
end