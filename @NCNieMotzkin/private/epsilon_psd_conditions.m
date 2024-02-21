function constraints = ...
    epsilon_psd_conditions(obj, state_real, state_img, mm, lm, epsilon)
%BASE_PSD_CONDITIONS Enforce essentially positive semi-definite moment matrices
% PARAMS:
%   state_real - The (real) state to constrain the matrices of.
%   state_img - The imaginary part of constraining state, or 'false'.
%      mm - The moment matrix (operators representation).
%      lm - The family of localizing matrices (operator representation).
% epsilon - Small positive number.

    if ~isequal(state_img, false)
        constraints = complex_epsilon_psd(obj, state_real, state_img, ...
                                          mm, lm, epsilon);
    else
        constraints = real_epsilon_psd(obj, state_real, mm, lm, epsilon);
    end
end

%% Private functions
function constraints = complex_epsilon_psd(obj, state_real, state_img, ...
                                           mm, lm, epsilon)

    espilon_mm = epsilon * eye(size(mm));
    espilon_lm = epsilon * eye(size(lm.max_sphere)); % [all lm same size]

    % PSD moment matrix
    constraints = [mm.Apply(state_real, state_img) + espilon_mm >= 0];

    % PSD localizing matrices
    if obj.exterior
        constraints = [constraints, ...
           lm.min_sphere.Apply(state_real, state_img) + espilon_lm >= 0];
    end
    constraints = [constraints, ...
       lm.max_sphere.Apply(state_real, state_img) + espilon_lm >= 0];

    for idx=1:3
        constraints = [constraints, ...
            lm.comm_plus{idx}.Apply(state_real, state_img) + espilon_lm >= 0];
        constraints = [constraints, ...
            lm.comm_minus{idx}.Apply(state_real, state_img) + espilon_lm >= 0];
    end
end

function constraints = real_epsilon_psd(obj, state_real, mm, lm, epsilon)

    espilon_mm = epsilon * eye(size(mm));
    espilon_lm = epsilon * eye(size(lm.max_sphere)); % [all lm same size]
    
    % PSD moment matrix
    constraints = [mm.Apply(state_real) + espilon_mm >= 0];

    % PSD localizing matrices
    if obj.exterior
        constraints = [constraints, ...
           lm.min_sphere.Apply(state_real) + espilon_lm >= 0];
    end
    constraints = [constraints, ...
       lm.max_sphere.Apply(state_real) + espilon_lm >= 0];

    for idx=1:3
        constraints = [constraints, ...
            lm.comm_plus{idx}.Apply(state_real) + espilon_lm >= 0];
        constraints = [constraints, ...
            lm.comm_minus{idx}.Apply(state_real) + espilon_lm >= 0];
    end
end
