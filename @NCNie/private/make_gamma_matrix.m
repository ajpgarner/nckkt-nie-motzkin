function gamma = make_gamma_matrix(obj, lm_level)
%MAKE_GAMMA_MATRIX Make Γ matrix for state optimality conditions.
% Γ = <f> - 1/2{f, MM}
%
    % Get LM of objective
    lm_f = obj.Scenario.LocalizingMatrix(obj.Objective, lm_level);

    % Get MM of appropriate level
    mm = obj.Scenario.MomentMatrix(lm_level);
    
    % Calculate anticommuting parts of matrix (manually)
    f_mm = (-0.5*obj.Objective) .* mm;
    mm_f = mm .* (-0.5*obj.Objective);
   
    % Sum together: 
    gamma = lm_f + f_mm + mm_f;
    
end

