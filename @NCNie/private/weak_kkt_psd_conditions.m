function c = weak_kkt_psd_conditions(obj, states, mm, lm)
%WEAK_KKT_PSD_CONDITIONS Create positivity constraints on matrices.
%   All PSD constraints are applied to sigma and mu
%

    % MM >= 0, all LM >= 0 for sigma
    c = base_psd_conditions(obj, states.sigma, false, mm, lm);
    
    % MM >= 0, all LM >= 0 for all mu
    for idx = 1:obj.d
        c = [c, base_psd_conditions(obj, states.mu.neg_square{idx}, false, mm, lm)];
    end
    for jdx = 1:(obj.d*obj.d)
        c = [c, base_psd_conditions(obj, states.mu.pos_comm{jdx}.a, ...
                                    states.mu.pos_comm{jdx}.b, mm, lm)];
        c = [c, base_psd_conditions(obj, states.mu.neg_comm{jdx}.a, ...
                                    states.mu.neg_comm{jdx}.b, mm, lm)];
    end
            
    % Lambda is almost unconstrained in this problem, because already 
    % (by architecture of Moment) lambda maps Hermitian operators to real
    % numbers, and lambda(a_i^2) = lambda(id), etc. since all operator
    % substitutions are applied before the moment is taken.    
    
end

