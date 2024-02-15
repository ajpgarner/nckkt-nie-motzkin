function kkt = mu_g_conditions(obj, mu)
%MU_G_CONDITIONS State constraints of the form μ_i(g_i) = 0

    kkt = [];
    
    % Constraints b_i^2 - 1
    for i = 1:obj.d 
        poly = obj.Constraints.neg_square(i);
        kkt = [kkt, poly.yalmip(mu.neg_square{i}) == 0];
    end
    
    % Constraints: δ ± i[a_j, b_k] (potentially complex states!)
    for idx = 1:(obj.d*obj.d)
        poly_pos = obj.Constraints.pos_comm(idx);
        poly_neg = obj.Constraints.neg_comm(idx);
        kkt = [kkt, poly_pos.yalmip(mu.pos_comm{idx}.a, mu.pos_comm{idx}.b) == 0, ...
                    poly_neg.yalmip(mu.neg_comm{idx}.a, mu.neg_comm{idx}.b) == 0];
    end
    
end

