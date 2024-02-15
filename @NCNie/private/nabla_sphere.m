function expr = nabla_sphere(obj, exterior, variate_idx, variate)
%NABLA_NEG_SQUARE  Derivative of sphere constraint.
% Interior constraint:
%   ∇(max - x_1^2 -x_2^2 -x_3^2) = -{x_1, y_1} - {x_2, y_2} - {x_3, y_3}.
% Exterior constriant:
%   ∇(x_1^2 + x_2^2 + x_3^2 - min) = {x_1, y_1} + {x_2, y_2} + {x_3, y_3}.
%
% PARAMS:
%     exterior - True for exterior of sphere (min_sphere constraint),
%                otherwise interior of sphere (max_sphere constraint).
%  variate_idx - The element of the variate that is non-zero (1, 2 or 3).
%      variate - The variate monomial.
%
   
   % Only one variate is non-zero, and it is x1, x2 or x3.
   switch variate_idx
       case 1
           expr = anticommutator(obj.x1, variate);
       case 2
           expr = anticommutator(obj.x2, variate);
       case 3
           expr = anticommutator(obj.x3, variate);
       otherwise
           error("Bad variate index.");
   end
   
   % Minus sign for interior constraint
   if ~exterior
       expr = -1.0 * expr;
   end
    
end

