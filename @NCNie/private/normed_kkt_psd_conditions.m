function c = normed_kkt_psd_conditions(obj, states, mm, lm)
%NORMED_KKT_PSD_CONDITIONS Create positivity constraints on matrices.
%   All PSD constraints are applied to all states.
%
   
    % MM >= 0, all LM >= 0 for sigma
    c = base_psd_conditions(obj, states.sigma, false, mm, lm);
    
    % MM >= 0, all LM >= 0 for all mu
    for idx = 1:obj.d
        c = [c, base_psd_conditions(obj, states.mu.neg_square{idx}, false, mm, lm)];
    end
    for jdx = 1:(obj.d*obj.d)
        c = [c, base_psd_conditions(obj, states.mu.pos_comm{jdx}.a, ...
                                         states.mu.pos_comm{jdx}.b, ...
                                         mm, lm)];
        c = [c, base_psd_conditions(obj, states.mu.neg_comm{jdx}.a, ...
                                         states.mu.neg_comm{jdx}.b, ...
                                         mm, lm)];
    end
        
    % MM >= 0, all LM >= 0 for all lambda+ and lambda-
    for kdx = 1:obj.d
        c = [c, base_psd_conditions(obj, states.lambda_plus{kdx}, false, mm, lm)];
        c = [c, base_psd_conditions(obj, states.lambda_minus{kdx}, false, mm, lm)];
    end        
end

