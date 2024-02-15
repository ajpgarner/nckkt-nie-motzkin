function [mm, lm] = make_matrices(obj, mm_level, lm_level)
%MAKE_MATRICES Make moment and localizing matrices

    % Make moment matrix
    mm = obj.Scenario.MomentMatrix(mm_level);
    
    % Make localizing matrices for 1 - x_1^2 - x_2^2 - x_3^2 >= 0
    lm = struct;
    lm.sphere = obj.Scenario.LocalizingMatrix(obj.Constraints.sphere, ...
                                              lm_level);
        
    % Make localizing matrices for  δ ± i[x_j, x_k] >= 0
    lm.comm_plus = cell(1, 3);   
    lm.comm_minus = cell(1, 3);
    for idx = 1:(obj.d*obj.d)
        lm.comm_plus{idx} = obj.Scenario.LocalizingMatrix(...
            obj.Constraints.comm_plus(idx), lm_level);
        lm.comm_minus{idx} = obj.Scenario.LocalizingMatrix(...
            obj.Constraints.comm_minus(idx), lm_level);                                      
    end
end

