function expr = nabla_comm(obj, a_idx, b_idx, pm, variate_idx, variate)
%NABLA_B  Derivative of expression: \delta \pm i[a_i, b_j]
%   =  \pm i [a_i, \bar{b}_j} + \mp i [b_j, \bar{a}_i]
%
% PARAMS:
%       a_idx - The subscript i of a_i.
%       b_idx - The subscript j of b_i.
%          pm - Prefactor in front of the commutator (set to +1 or -1).
% variate_idx - The element of the variate that is non-zero (from 1 to 2d).
%     variate - The variate monomial.
%
   
   % Only one variate term is non-zero:
   if variate_idx == a_idx % \mp i [b_i, \bar[b]_i]
       expr = (-1i * pm) * commutator(obj.Scenario.get(obj.d+b_idx), variate);
   elseif variate_idx == obj.d + b_idx % \pm i [a_i, \bar[b]_j]
       expr = (1i * pm) * commutator(obj.Scenario.get(a_idx), variate);
   else % Otherwise, zero:
       expr = obj.Scenario.zero;
   end
    
end