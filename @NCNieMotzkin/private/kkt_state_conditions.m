function kkt = kkt_state_conditions(obj, states, kkt_level)
%KKT_STATE_CONDITIONS Scalar conditions for normed KKT.

    % mu_i(g_i) == 0 state conditions (here, mostly mu(sg) = mu(gs) == 0).
    kkt = mu_g_conditions(obj, states.mu, kkt_level);

    % Optimality conditions for each polynomial:
    monomials = obj.Scenario.WordList(kkt_level);
    for vIdx = 1:3 % Y1, Y2, Y3:
        for mIdx = 1:numel(monomials)
            monomial = monomials(mIdx);
            kkt = [kkt, optimality_condition(obj, states, vIdx, monomial)];
        end
    end
end

%% Private functions
function kkt = mu_g_conditions(obj, mu, kkt_level)
%MU_G_CONDITIONS State constraints of the form μ_i(g_i) = 0
%

    % Since we are doing essential conditions, we need our 's'
    monomials = obj.Scenario.WordList(kkt_level);

    %  Constraint: max - x_1^2 - x_2^2 - x_3^2 >= 0
    kkt = essential_mu_g_s(obj.Constraints.max_sphere,...
                           monomials, mu.max_sphere);
        
    %  Constraint: x_1^2 + x_2^2 + x_3^2 - min >= 0
    if obj.exterior
    % TODO: Special essential case of min_sphere    
        kkt = [kkt, essential_mu_g_s(obj.Constraints.min_sphere, ...
                                     monomials, mu.min_sphere)];
    end
       
    % Constraints: δ ± i[x_j, x_k] >= 0
    for idx = 1:3
        kkt = [kkt, essential_mu_g_s(obj.Constraints.comm_plus{idx}, ...
                                     monomials, mu.comm_plus{idx}), ...
                    essential_mu_g_s(obj.Constraints.comm_minus{idx}, ...
                                     monomials, mu.comm_minus{idx})];
    end
end

function result = essential_mu_g_s(constraint, monomials, state)
        
    % Test if state is complex:
    complex = isstruct(state);
    
    % g * s (broadcast; filter to elements that exist in MM/LMs)
    gs = constraint .* monomials;
    gs = gs.onlyExistingSymbols;
    if ~isempty(gs)
        if complex
            result = [gs.yalmip(state.a, state.b) == 0];
        else
            result = [gs.yalmip(state) == 0];
        end
    else
        result = [];
    end
    
    % s *g (broadcast; filter to elements that exist in MM/LMs)
    sg = monomials .* constraint;
    sg = sg.onlyExistingSymbols;
    if ~isempty(sg)
        if complex
             result = [result, sg.yalmip(state.a, state.b) == 0];
        else
             result = [result, sg.yalmip(state) == 0];
        end
    end    
end
