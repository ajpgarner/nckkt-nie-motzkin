function condition = optimality_condition(obj, states, vIdx, variate)
%OPTIMALITY_CONDITION Generates a single optimality condition.
%
% Since there are no equality constraints, this would be the same whether
% normed, weak or essential.
%
% PARAMS:
%    states - The structure of yalmip sdpvars
%    vIdx   - The index of the non-zero term in the variate vector
%  variate  - The variate monomial to be substituted in for the vIdx term.
%

    % Print name of variate in verbose mode
    if obj.Verbose >=2
        fprintf("\ny_%d = %s:\n", vIdx, variate.ObjectName);        
    end
    
    % Objective function
    nabla_f = obj.nabla_objective(vIdx, variate);
    expr = nabla_f.Apply(states.sigma);    
    if obj.Verbose >= 2
        fprintf(" ∇f = %s\n", nabla_f.ObjectName);
    end
 
    % Inequality constraint: Sphere maximum
    D_max_sphere = nabla_sphere(obj, vIdx, variate);
    if ~D_max_sphere.IsZero
        expr = expr - D_max_sphere.Apply(states.mu.max_sphere);
    end
    if obj.Verbose>=2 
        fprintf(" ∇max-sphere = %s\n", D_max_sphere.ObjectName);
    end
    
    % Inequality constraint: Sphere minimum (in exterior mode)
    % NB: D_min_sphere = - D_max_sphere
    if obj.exterior && ~D_max_sphere.IsZero
        expr = expr + D_max_sphere.Apply(states.mu.min_sphere);
        if obj.Verbose>=2 
            D_min_sphere = -D_max_sphere;
            fprintf(" ∇min-sphere = %s\n", D_min_sphere.ObjectName);
        end
    end
       
    % Inequality constraints: Commutator expressions
    expr = commutator_constraints(obj, expr, vIdx, variate, states.mu);
    
    % Convert expression to an equality:
    condition = [expr == 0];
end

%% Private functions
function expr = commutator_constraints(obj, expr, vIdx, variate, mu)
   % In principle there are six constraints, but most terms will be zero:
   %  #1: ∇(δ ± i[x2, x3]) = ± i[x2, y3] ∓ i[x3, y2]
   %  #2: ∇(δ ± i[x3, x1]) = ± i[x3, y1] ∓ i[x1, y3]
   %  #3: ∇(δ ± i[x1, x2]) = ± i[x1, y2] ∓ i[x2, y1]
   
   % Non-zero terms depend on variate index:
   switch vIdx
       case 1
           % #2: ∇(δ ± i[x3, x1]) = ± i[x3, y1] ∓ i[x1, y3], first term
           comm_first  = 1i * commutator(obj.x3, variate);
           pIdx = 2;
           
           % #3: ∇(δ ± i[x1, x2]) = ± i[x1, y2] ∓ i[x2, y1], second term
           comm_second = -1i * commutator(obj.x2, variate);
           mIdx = 3;
           
       case 2
           % #3: ∇(δ ± i[x1, x2]) = ± i[x1, y2] ∓ i[x2, y1], first term
           comm_first = 1i * commutator(obj.x1, variate);
           pIdx = 3;
           
           % #1: ∇(δ ± i[x2, x3]) = ± i[x2, y3] ∓ i[x3, y2], second term
           comm_second = -1i * commutator(obj.x3, variate);
           mIdx = 1;
           
       case 3
           % #1: ∇(δ ± i[x2, x3]) = ± i[x2, y3] ∓ i[x3, y2], first term
           comm_first = 1i * commutator(obj.x2, variate);
           pIdx = 1;
           
           % #2: ∇(δ ± i[x3, x1]) = ± i[x3, y1] ∓ i[x1, y3], second term
           comm_second = -1i * commutator(obj.x1, variate);
           mIdx = 2;           
   end
   
   % Add first terms if non-trivial
   if ~comm_first.IsZero
       expr = expr - comm_first.Apply(mu.comm_plus{pIdx}.a, ...
                                       mu.comm_plus{pIdx}.b);
       expr = expr + comm_first.Apply(mu.comm_minus{pIdx}.a, ...
                                       mu.comm_minus{pIdx}.b);    
   end
   
   % Add second terms if non-tirival
   if ~comm_second.IsZero       
       expr = expr - comm_second.Apply(mu.comm_plus{mIdx}.a, ...
                                        mu.comm_plus{mIdx}.b);
       expr = expr + comm_second.Apply(mu.comm_minus{mIdx}.a, ...
                                        mu.comm_minus{mIdx}.b);
   end
   
   % Output derivatives in verbose mode:
   if obj.Verbose >= 2
       if ~comm_first.IsZero
           fprintf(" ∇(δ+i[x%d,x%d]) = %s\n", ...
                   pIdx, vIdx, comm_first.ObjectName);
           neg_first = -comm_first;
           fprintf(" ∇(δ-i[x%d,x%d]) = %s\n", ...
                   pIdx, vIdx, neg_first.ObjectName);
       end
       if ~comm_second.IsZero
           fprintf(" ∇(δ+i[x%d,x%d]) = %s\n", ...
                   vIdx, mIdx, comm_second.ObjectName);
           neg_second = -comm_second;
           fprintf(" ∇(δ-i[x%d,x%d]) = %s\n", ...
                   vIdx, mIdx, neg_second.ObjectName);
       end
   end
end