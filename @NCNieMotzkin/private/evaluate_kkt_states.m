function val_states = evaluate_kkt_states(obj, ym_states)
%EVALUATE_NORMED_KKT_STATES Get actual numbers associated with solution.

    val_states = struct;
    val_states.sigma = value(ym_states.sigma);
    
    % Collection of mu states
    val_states.mu = evaluate_mu(obj, ym_states.mu);
    
end

%% Private functions
function mu = evaluate_mu(obj, ym_states)
    mu = struct;
    if obj.exterior
        mu.min_sphere = value(ym_states.min_sphere);
    end
    
    mu.max_sphere = value(ym_states.max_sphere);
    
    mu.comm_plus = evaluate_complex_states(3, ym_states.comm_plus);
    mu.comm_minus = evaluate_complex_states(3, ym_states.comm_minus);
end

function res = evaluate_complex_states(N, ym_states)
    res = cell(1, N);
    for i = 1:N
        res{i} = struct;
        res{i}.a = value(ym_states{i}.a);
        res{i}.b = value(ym_states{i}.b);
    end
end