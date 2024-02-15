function kkt = kkt_state_conditions(obj, states, kkt_level)
%KKT_STATE_CONDITIONS Scalar conditions for normed KKT.

    % mu_i(g_i) = 0 state conditions
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

    %  Constraint: max - x_1^2 - x_2^2 - x_3^2 >= 0
    kkt = [obj.Constraints.max_sphere.yalmip(mu.max_sphere.a, mu.max_sphere.b) == 0];
    
    %  Constraint: x_1^2 + x_2^2 + x_3^2 - min >= 0
    if obj.exterior
    % TODO: Special essential case of min_sphere    
        kkt = [kkt, essential_mu_g_s_conditions(obj, mu.min_sphere, kkt_level)];
    end
       
    % Constraints: δ ± i[x_j, x_k] >= 0
    for idx = 1:3
        kkt = [kkt, obj.Constraints.comm_plus{idx}.yalmip(...
                        mu.comm_plus{idx}.a, mu.comm_plus{idx}.b) == 0, ...
                    obj.Constraints.comm_minus{idx}.yalmip(...
                        mu.comm_minus{idx}.a, mu.comm_minus{idx}.b) == 0];
    end
end

function result = essential_mu_g_s_conditions(obj, state, kkt_level)   
    monomials = obj.Scenario.WordList(kkt_level);
    
    % g * s (broadcast; filter to elements that exist in MM/LMs)
    gs = obj.Constraints.min_sphere .* monomials;
    gs = gs.onlyExistingSymbols;
    if ~isempty(gs)
        result = [gs.yalmip(state.a, state.b) == 0];
    else
        result = [];
    end
    
    % s *g (broadcast; filter to elements that exist in MM/LMs)
    sg = monomials .* obj.Constraints.min_sphere;
    sg = sg.onlyExistingSymbols;
    if ~isempty(sg)
        result = [result, sg.yalmip(state.a, state.b) == 0];
    end    
end
