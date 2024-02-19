function constraints = base_psd_conditions(obj, state_real, state_img, mm, lm)
%BASE_PSD_CONDITIONS Enforce positive semi-definite moment matrices
% PARAMS:
%   state_real - The (real) state to constrain the matrices of.
%   state_img - The imaginary part of constraining state, or 'false'.
%      mm - The moment matrix (operators representation).
%      lm - The family of localizing matrices (operator representation).
%

    if ~isequal(state_img, false) % Do we have imaginary part too?
        constraints = complex_psd(obj, state_real, state_img, mm, lm);
    else
        constraints = real_psd(obj, state_real, mm, lm);
    end    
end



%% Private functions
function constraints = complex_psd(obj, state_real, state_img, mm, lm)

    % PSD moment matrix
    constraints = [mm.yalmip(state_real, state_img) >= 0];

    % PSD localizing matrices
    if obj.exterior
        constraints = [constraints, ...
           lm.min_sphere.yalmip(state_real, state_img) >= 0];
    end
    constraints = [constraints, ...
       lm.max_sphere.yalmip(state_real, state_img) >= 0];

    for idx=1:3
        constraints = [constraints, ...
            lm.comm_plus{idx}.yalmip(state_real, state_img) >= 0];
        constraints = [constraints, ...
            lm.comm_minus{idx}.yalmip(state_real, state_img) >= 0];
    end
end

function constraints = real_psd(obj, state_real, mm, lm)
    % PSD moment matrix
    constraints = [mm.yalmip(state_real) >= 0];

    % PSD localizing matrices
    if obj.exterior
        constraints = [constraints, lm.min_sphere.yalmip(state_real) >= 0];
    end
    constraints = [constraints, lm.max_sphere.yalmip(state_real) >= 0];

    for idx=1:3
        constraints = [constraints, ...
                       lm.comm_plus{idx}.yalmip(state_real) >= 0];
        constraints = [constraints, ...
                       lm.comm_minus{idx}.yalmip(state_real) >= 0];
    end
end
