function [mm_eps, lm_eps] = make_epsilon_matrices(obj, mm, lm, epsilon)
%MAKE_EPSILON_MATRICES Add epsilon to MM and LM

    % Validate epsilon
    assert(isscalar(epsilon) && (epsilon > 0), ...
           "Epsilon must be a positive number.");
    
    % Make epsilon perturbation matrix
    mm_dim = size(mm, 1);
    mm_eye_eps = MTKValueMatrix(obj.Scenario, epsilon * eye(mm_dim));
    mm_eps = mm + mm_eye_eps;
    
    
    % Make localizing matrices for 1-b_i^2 >= 0
    lm_eps = struct;
    lm_eps.neg_square = cell(1, obj.d);
    lm_eps.pos_comm = cell(1, obj.d * obj.d);   
    lm_eps.neg_comm = cell(1, obj.d * obj.d);
    
    % Also create LM ID matrix if different size from MM ID matrix
    lm_dim = size(lm.neg_square{1}, 1);
    if lm_dim ~= mm_dim
        lm_eye_eps = MTKValueMatrix(obj.Scenario, epsilon * eye(lm_dim));        
    else
        lm_eye_eps = mm_eye_eps;
    end
    
    % Make espilon perturned LM
    for idx = 1:obj.d
        lm_eps.neg_square{idx} = lm.neg_square{idx} + lm_eye_eps;        
    end
    for jdx = 1:(obj.d*obj.d)        
        lm_eps.pos_comm{jdx} = lm.pos_comm{jdx} + lm_eye_eps;
        lm_eps.neg_comm{jdx} = lm.neg_comm{jdx} + lm_eye_eps;
    end   
    
end

