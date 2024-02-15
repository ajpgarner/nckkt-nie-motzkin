function kkt = weak_kkt_state_conditions(obj, states, kkt_level)
%WEAK_KKT_STATE_CONDITIONS Scalar conditions for weak KKT.
    
    % mu_i(g_i) = 0 state conditions
    kkt = mu_g_conditions(obj, states.mu);
        
    % Optimality conditions for each polynomial:
    monomials = obj.Scenario.WordList(kkt_level);
    for vIdx = 1:(obj.d*2)
        for mIdx = 1:numel(monomials)
            monomial = monomials(mIdx);
            kkt = [kkt, ...
                   weak_kkt_optimality_condition(obj, states, ...
                                                   vIdx, monomial)];
        end        
    end   
end