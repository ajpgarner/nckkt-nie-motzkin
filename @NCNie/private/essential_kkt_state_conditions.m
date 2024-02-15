function kkt = essential_kkt_state_conditions(obj, states, kkt_level)
%ESSENTIAL_KKT_STATE_CONDITIONS Scalar conditions for weak KKT.
    
    % mu_i(g_i) = 0 state conditions
    kkt = mu_g_conditions(obj, states.mu);
    
    % Calculate word list
    monomials = obj.Scenario.WordList(kkt_level);
        
    % Additional essential mu_i(g_i s) = 0, mu_i(s g_i) = 0 conditions:
    kkt = [kkt, mu_g_s_conditions(obj, states, monomials)];
            
    % Optimality conditions for each polynomial:    
    for vIdx = 1:(obj.d*2)
        for mIdx = 1:numel(monomials)
            monomial = monomials(mIdx);
            % Essential optimality cond. same as weak optimality cond.
            kkt = [kkt, weak_kkt_optimality_condition(obj, states, ...
                                                      vIdx, monomial)];
        end        
    end   
end

%% Private functions
function c = mu_g_s_conditions(obj, states, monomials)
    c = [];
    for idx=1:obj.d
        c = [c, add_gs_sg_constraint(states.mu.neg_square{idx}, ...
                                     monomials, ...
                                     obj.Constraints.neg_square(idx))];        
    end

    for jdx=1:(obj.d*obj.d)
        c = [c, add_complex_gs_sg_constraint(states.mu.pos_comm{jdx}.a, ...
                                     states.mu.pos_comm{jdx}.b, ...
                                     monomials, ...
                                     obj.Constraints.pos_comm(jdx))];
        c = [c, add_complex_gs_sg_constraint(states.mu.neg_comm{jdx}.a, ...
                                     states.mu.neg_comm{jdx}.b, ...
                                     monomials, ...
                                     obj.Constraints.neg_comm(jdx))];                                         
    end
end

function result = add_gs_sg_constraint(state, monomials, constraint)   
    % g * s (broadcast; filter to elements that exist in MM/LMs)
    gs = constraint .* monomials;
    gs = gs.onlyExistingSymbols;
    if ~isempty(gs)
        result = [gs.yalmip(state) == 0];
    else
        result = [];
    end
    
    % s *g (broadcast; filter to elements that exist in MM/LMs)
    sg = monomials .* constraint;
    sg = sg.onlyExistingSymbols;
    if ~isempty(sg)
        result = [result, sg.yalmip(state) == 0];
    end    
end


function result = add_complex_gs_sg_constraint(state_real, state_img, monomials, constraint)   
    % g * s (broadcast; filter to elements that exist in MM/LMs)
    gs = constraint .* monomials;
    gs = gs.onlyExistingSymbols;
    if ~isempty(gs)
        result = [gs.yalmip(state_real, state_img) == 0];
    else
        result = [];
    end
    
    % s *g (broadcast; filter to elements that exist in MM/LMs)
    sg = monomials .* constraint;
    sg = sg.onlyExistingSymbols;
    if ~isempty(sg)
        result = [result, sg.yalmip(state_real, state_img) == 0];
    end    
end