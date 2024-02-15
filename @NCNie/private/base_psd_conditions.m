function constraints = base_psd_conditions(obj, state_real, state_img, mm, lm)
%BASE_PSD_CONDITIONS Enforce positive semi-definite moment matrices
% PARAMS:
%   state_real - The (real) state to constrain the matrices of.
%   state_img - The imaginary part of constraining state, or 'false'.
%      mm - The moment matrix (operators representation).
%      lm - The family of localizing matrices (operator representation).
%

    if ~isequal(state_img, false) % Do we have imaginary part too?
        % PSD moment matrix
        constraints = [mm.yalmip(state_real, state_img) >= 0];

        % PSD localizing matrices
        for idx=1:(obj.d)
            constraints = [constraints, ...
                lm.neg_square{idx}.yalmip(state_real, state_img) >= 0];
        end
        for idx=1:(obj.d*obj.d)
            constraints = [constraints, ...
                lm.pos_comm{idx}.yalmip(state_real, state_img) >= 0];
            constraints = [constraints, ...
                lm.neg_comm{idx}.yalmip(state_real, state_img) >= 0];
        end
        
    else % Purely real:        
        % PSD moment matrix
        constraints = [mm.yalmip(state_real) >= 0];

        % PSD localizing matrices
        for idx=1:(obj.d)
            constraints = [constraints, ...
                lm.neg_square{idx}.yalmip(state_real) >= 0];
        end
        for idx=1:(obj.d*obj.d)
            constraints = [constraints, ...
                lm.pos_comm{idx}.yalmip(state_real) >= 0];
            constraints = [constraints, ...
                lm.neg_comm{idx}.yalmip(state_real) >= 0];
        end
    end
end
