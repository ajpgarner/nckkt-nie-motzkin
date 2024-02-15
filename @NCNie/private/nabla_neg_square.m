function expr = nabla_neg_square(obj, op_idx, variate_idx, variate)
%NABLA_NEG_SQUARE  Derivative of expression: 1 - x_i^2 = -{x_i, \bar{x}_i},
% where x = [a_1, ... a_d, b_1, ... b_d]
%
% PARAMS:
%       op_idx - The operator to square
%  variate_idx - The element of the variate that is non-zero (from 1 to 2d).
%      variate - The variate monomial.
%
   
   % Only one variate is non-zero:
   if variate_idx == op_idx
       expr = -1.0 * anticommutator(obj.Scenario.get(op_idx), variate);
   else
       expr = obj.Scenario.zero;
   end
    
end

