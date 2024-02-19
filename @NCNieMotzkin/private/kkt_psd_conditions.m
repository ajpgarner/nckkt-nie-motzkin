function c = kkt_psd_conditions(obj, states, mm, lm, epsilon)
% KKT_PSD_CONDITIONS Create positivity constraints on matrices for all
% sigma and mu.
%
% NB: In exterior mode, essential positivity is applied for min_sphere.
%

    % MM >= 0, all LM >= 0 for sigma
    c = base_psd_conditions(obj, states.sigma, false, mm, lm);
    
    % MM >= 0, all LM >= 0 for interior sphere mu, with essential epsilon.
    c = [c, epsilon_psd_conditions(obj, states.mu.max_sphere, false, ...
                                   mm, lm, epsilon)];

    % If min sphere, then add essential PSD constraint for associated mu
    if obj.exterior
        c = [c, epsilon_psd_conditions(obj, states.mu.min_sphere, false,...
                                       mm, lm, epsilon)];
    end
                                    
    % MM >= 0, all LM >= 0 for each commutator mu, with essential epsilon.
    for idx = 1:3
        c = [c, epsilon_psd_conditions(obj, states.mu.comm_plus{idx}.a, ...
                                            states.mu.comm_plus{idx}.b, ...
                                            mm, lm, epsilon)];
        c = [c, epsilon_psd_conditions(obj, states.mu.comm_minus{idx}.a,...
                                            states.mu.comm_minus{idx}.b,...
                                            mm, lm, epsilon)];
    end
end
