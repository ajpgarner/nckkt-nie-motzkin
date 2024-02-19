function expr = nabla_sphere(obj, variate_idx, variate)
%NABLA_SPHERE Derivative of max sphere constraint.
% Interior (max sphere) constraint:
%   ∇(max - x_1^2 -x_2^2 -x_3^2) = -{x_1, y_1} - {x_2, y_2} - {x_3, y_3}.
% Exterior (min sphere) constriant:
%   ∇(x_1^2 + x_2^2 + x_3^2 - min) = {x_1, y_1} + {x_2, y_2} + {x_3, y_3}.
%
% Both derivatives are independent of value of min/max and differ by sign.
% This function returns max sphere; for min sphere, multiply by -1.
%
% PARAMS:
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
       
end

