function condition = normed_kkt_optimality_condition(obj, states, vIdx, variate)
%NORMED_KKT_OPTIMALITY_CONDITION Generates normed optimality condition.
%
% PARAMS:
%    states - The structure of yalmip sdpvars
%    vIdx   - The index of the non-zero term in the variate vector
%  variate  - The variate monomial to be substituted in for the vIdx term.
%

    %% Print name of variate in verbose mode
    if obj.Verbose >=2
        if vIdx <= obj.d
            v_prefix = "a";
            v_num = vIdx;
        else
            v_prefix = "b";
            v_num = vIdx - obj.d;
        end
        header_str = sprintf("\\bar{%s}_%d = %s", ...
                v_prefix, v_num, variate.ObjectName);
        
        fprintf("\n");
        fprintf("%s:\n", header_str);        
    end
    
    %% Objective
    nabla_f = obj.nabla_objective(vIdx, variate);
    if obj.Verbose >= 2
        fprintf(" ∇f = %s\n", nabla_f.ObjectName);
    end

    expr = nabla_f.yalmip(states.sigma);    
 
    %% Inequality constraint:  b_^2 - 1 >= 0
    if vIdx > obj.d % Only has support on in \bar{b}_i term
        D_neg_bb = obj.nabla_neg_square(vIdx, vIdx, variate);
        if ~D_neg_bb.IsZero
            expr = expr - D_neg_bb.yalmip(states.mu.neg_square{vIdx-obj.d});            
        end
        if obj.Verbose>=2 
            fprintf(" ∇(1-b_%d^2) = %s\n", ...
                    (vIdx-obj.d), D_neg_bb.ObjectName);
        end
    end
    
    %% Inequality constraints: Commutator expressions
    if vIdx <= obj.d  % \bar{a}_i variate: loop over j for [a_i, b_j]
        for j = 1:obj.d
            joint_idx = (vIdx-1)*obj.d + j; % 'i*d + j' (up to 1-indexing)
            
            D_pos_comm = obj.nabla_comm(vIdx, j, +1, vIdx, variate);
            if ~D_pos_comm.IsZero
                expr = expr ...
                     - D_pos_comm.yalmip(states.mu.pos_comm{joint_idx}.a, ...
                                         states.mu.pos_comm{joint_idx}.b);
            end
                
            D_neg_comm = obj.nabla_comm(vIdx, j, -1, vIdx, variate);
            if ~D_neg_comm.IsZero
                expr = expr ...
                     - D_neg_comm.yalmip(states.mu.neg_comm{joint_idx}.a, ...
                                         states.mu.neg_comm{joint_idx}.b);
            end
            
            if obj.Verbose>=2 
                fprintf(" ∇(+i[a%d, b%d]) = %s\n", ...
                        vIdx, j, D_pos_comm.ObjectName);
                fprintf(" ∇(-i[a%d, b%d]) = %s\n", ...
                        vIdx, j, D_pos_comm.ObjectName);
            end
        end
    else % \bar{b}_j variate: loop over i for [a_i, b_j]
        for i = 1:obj.d
            joint_idx = (i-1)*obj.d + (vIdx - obj.d);
            
            D_pos_comm = obj.nabla_comm(i, vIdx - obj.d, +1, vIdx, variate);
            if ~D_pos_comm.IsZero
                expr = expr ...
                     - D_pos_comm.yalmip(states.mu.pos_comm{joint_idx}.a, ...
                                         states.mu.pos_comm{joint_idx}.b);
            end
                        
            D_neg_comm = obj.nabla_comm(i, vIdx - obj.d, -1, vIdx, variate);
            if ~D_neg_comm.IsZero
                expr = expr ...
                     - D_neg_comm.yalmip(states.mu.neg_comm{joint_idx}.a, ...
                                         states.mu.neg_comm{joint_idx}.b);
            end
            
           if obj.Verbose>=2 
                fprintf(" ∇(+i[a%d, b%d]) = %s\n", ...
                        i, vIdx - obj.d, D_pos_comm.ObjectName);
                fprintf(" ∇(-i[a%d, b%d]) = %s\n", ...
                        i, vIdx - obj.d, D_pos_comm.ObjectName);
            end
        end
    end
    
    
    %% Equality constraints: a_i^2 - 1 = 0    
    if vIdx <= obj.d % Only has support on \bar{a}_i term
        D_neg_aa = obj.nabla_neg_square(vIdx, vIdx, variate);
        if ~D_neg_aa.IsZero
            expr = expr - D_neg_aa.yalmip(states.lambda_plus{vIdx}) ...
                        + D_neg_aa.yalmip(states.lambda_minus{vIdx});       
        end
        if obj.Verbose>=2 
            fprintf(" ∇(1-a_%d^2) = %s\n", ...
                    vIdx, D_neg_aa.ObjectName);
        end
    end
  
    % Convert expression to an equality:
    condition = [expr == 0];
end