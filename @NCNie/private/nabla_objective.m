function expr = nabla_objective(obj, variate_idx, variate)
%NABLA_NEG_SQUARE  Derivative of objective function.
%
% We will use linearly:
%    ∇a_i = \bar{a}_i
%    ∇b_j = \bar{b}_j
%    ∇{a_i, b_j} = {a_i, \bar{b}_j} + {b_j, \bar{a}_i}
% PARAMS:
%  variate_idx - The element of the variate that is non-zero (from 1 to 2d).
%      variate - The variate monomial.
%

    FC = obj.ObjectiveTensor;

    % Does variate match a or b? Guaranteed to only match one...
    if variate_idx <= obj.d %  Case: \bar{a}_i

        % First, linear terms in objective:    
        if FC(variate_idx + 1, 1) ~= 0
            % Nabla k a_i = k \bar{a}_i
            expr = FC(variate_idx + 1, 1) * variate;
        else
            expr = obj.Scenario.zero;
        end

        % Since variate is substituted into some a_i,
        %  we need the terms {b_j, \bar{a}_i} for all j.
        for j = 1:obj.d
            k = 0.5 * FC(variate_idx + 1, j + 1);
            if k ~= 0
                expr = expr + k * anticommutator(obj.b(j), variate);
            end
        end
    else % Case : \bar{b}_j

        % First, linear terms in objective:    
        if FC(1, variate_idx - obj.d + 1) ~= 0
            % Nabla k b_j = k \bar{b}_j
            expr = FC(1, variate_idx - obj.d + 1) * variate;
        else
            expr = obj.Scenario.zero;
        end

        % Since variate is substituted into some b_j,
        %  we need the terms {a_i, \bar{b}_j} for all i.
        for i = 1:obj.d
            k = 0.5 * FC(i + 1, variate_idx - obj.d + 1);
            if k ~= 0
                expr = expr + k * anticommutator(obj.a(i), variate);
            end
        end
    end
end