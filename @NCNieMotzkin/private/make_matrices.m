function [mm, lm, gamma] = make_matrices(obj, mm_level, lm_level)
%MAKE_MATRICES Make moment and localizing matrices

    % Make gamma matrix?
    make_gamma = nargout >= 3;
    
    % Make moment matrix
    mm = obj.Scenario.MomentMatrix(mm_level);
    
    % Make localizing matrices for 1 - x_1^2 - x_2^2 - x_3^2 >= 0
    lm = struct;
    if obj.exterior
        lm.min_sphere = obj.Scenario.LocalizingMatrix(...
            obj.Constraints.min_sphere, lm_level);
    end
    lm.max_sphere = obj.Scenario.LocalizingMatrix(...
        obj.Constraints.max_sphere, lm_level);
        
    % Make localizing matrices for  δ ± i[x_j, x_k] >= 0
    lm.comm_plus = cell(1, 3);   
    lm.comm_minus = cell(1, 3);
    for idx = 1:3
        lm.comm_plus{idx} = obj.Scenario.LocalizingMatrix(...
            obj.Constraints.comm_plus{idx}, lm_level);
        lm.comm_minus{idx} = obj.Scenario.LocalizingMatrix(...
            obj.Constraints.comm_minus{idx}, lm_level);                                      
    end
    
    % Make gamma matrix (for state optimality purposes), if requested...
    if make_gamma
        gamma = make_gamma_matrix(obj, lm_level);
    end
end

%% Private functions
function gamma = make_gamma_matrix(obj, lm_level)
%MAKE_GAMMA_MATRIX Make Γ matrix for state optimality conditions.
% Γ = <f> - 1/2{f, MM}
%
    % Get LM of objective
    lm_f = obj.Scenario.LocalizingMatrix(obj.Objective, lm_level);

    % Get MM at same level as LM
    mm = obj.Scenario.MomentMatrix(lm_level);
    
    % Calculate anticommuting parts of matrix (manually)
    minus_half_obj = -0.5*obj.Objective;
    f_mm = minus_half_obj .* mm;
    mm_f = mm .* minus_half_obj;
   
    % Sum together: 
    gamma = lm_f + f_mm + mm_f;
    
end


