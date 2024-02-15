function conditions = ...
    state_optimality_conditions(obj, gamma, sigma, so_level)
%STATE_OPTIMALITY_CONDITIONS
%

    % Get commutating constraints...
    monomials = obj.Scenario.WordList(so_level);
    commutators = 1i * commutator(monomials, obj.Objective);
    commutators = commutators.onlyExistingSymbols; % Ignore non-existing expressions.
    
    % Make yalmipified conditions
    conditions = [commutators.yalmip(sigma) == 0];
       
    % Make PSD gamma-matrix condition
    ym_gamma = gamma.yalmip(sigma);
    ym_gamma = ym_gamma(2:end, 2:end);   
    conditions = [conditions, ym_gamma  >= 0];
    
end

