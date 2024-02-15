function c = ...
    essential_kkt_psd_conditions(obj, states, mm, lm, mm_eps, lm_eps)
%ESSENTIAL_KKT_PSD_CONDITIONS Create positivity constraints on matrices.
%   All PSD constraints are applied to sigma and mu.
% For essential case, on mu, the MM and LM have a small positive offset.
%

    % MM >= 0, all LM >= 0 for sigma
    c = base_psd_conditions(obj, states.sigma, false, mm, lm);
    
    % MM + eI >= 0, all LM + eI >= 0 for all mu
    for idx = 1:obj.d
        c = [c, base_psd_conditions(obj, states.mu.neg_square{idx}, false, ...
                                    mm_eps, lm_eps)];
    end
    for jdx = 1:(obj.d*obj.d)
        c = [c, base_psd_conditions(obj, states.mu.pos_comm{jdx}.a, ...
                                    states.mu.pos_comm{jdx}.b, ...
                                    mm_eps, lm_eps)];
        c = [c, base_psd_conditions(obj, states.mu.neg_comm{jdx}.a, ...
                                    states.mu.neg_comm{jdx}.b, ...
                                    mm_eps, lm_eps)];
    end
        
    % Lambda is almost unconstrained in this problem, because already 
    % (by architecture of Moment) lambda maps Hermitian operators to real
    % numbers, and lambda(a_i^2) = lambda(id), etc. since all operator
    % substitutions are applied before the moment is taken.    
    
end

