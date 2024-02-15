function expr = nabla_objective(obj, vIdx, variate)
%NABLA_OBJECTIVE Derivative of objective function.
%
% Linear components (writing variate vector as y_1, y_2, y_3):
%  ∇({x_1^4, x_2^2}) = {x_2^2, {x_1^2, {x_1, y_1}}} + {x_1^4, {x_2, y_2}}
%  ∇({x_1^2, x_2^4}) = {x_1^2, {x_2^2, {x_2, y_2}}} + {x_2^4, {x_1, y_1}}
%  ∇(x_3^6) = {x_3^4, {x_3, y_3}} + x_3^2 {x_3, y_3} x_3^2
%  ∇(x_1^2 x_2^2 x_3^2) = x_1^2 x_2^2 {x_3, y_3} + x_1^2 {x_2,y_2} x_3^2 
%                         + {x_1, y_1} x_2^2 x_3^2 
%  ∇(x_3^2 x_2^2 x_1^2) = x_3^2 x_2^2 {x_1, y_1} + x_3^2 {x_2,y_2} x_1^2 
%                         + {x_3, y_3} x_2^2 x_1^2 
%
% (NB: Objective also has real scaling factors!)
%
% PARAMS:
%     vIdx - The element of the variate that is non-zero (1, 2 or 3).
%  variate - The variate monomial.
%


    % Specialize behaviour depending on the variate index
    % (as many terms will be zero):
    switch vIdx
        case 1
            expr = nabla_objective_y1(obj, variate);         
        case 2
            expr = nabla_objective_y2(obj, variate);            
        case 3
            expr = nabla_objective_y3(obj, variate);
        otherwise
            error("Bad vIdx.");            
    end    
end

%% Private functions
function expr = nabla_objective_y1(obj, variate)
%NABLA_OBJECTIVE_Y1 y_1 specialization of derivative.

    % {x_1, y_1}:
    ac_x1_y1 = anticommutator(obj.x1, variate);  
    
    % First term of ∇({x_1^4, x_2^2}) -> 0.5 * {x_2^2, {x_1^2, {x_1, y_1}}}
    expr = 0.5 * anticommutator(obj.x2_pow2, ...
                                anticommutator(obj.x1_pow2, ac_x1_y1));
    
    % Second term of  ∇({x_1^2, x_2^4}) -> 0.5 * {x_2^4, {x_1, y_1}}
    expr = expr + 0.5 * anticommutator(obj.x2_pow4, ac_x1_y1);
    
    % Third term of ∇(x_1^2 x_2^2 x_3^2) -> -1.5 * {x_1, y_1} x_2^2 x_3^2 
    expr = expr - 1.5 * ac_x1_y1 * obj.x2_pow2 * obj.x3_pow2;
    
    % First term of ∇(x_3^2 x_2^2 x_1^2) -> -1.5 * x_3^2 x_2^2 {x_1, y_1}
    expr = expr - 1.5 * obj.x3_pow2 * obj.x2_pow2 * ac_x1_y1;
end

function expr = nabla_objective_y2(obj, variate) 
%NABLA_OBJECTIVE_Y2 y_2 specialization of derivative.

    % {x_2, y_2}:
    ac_x2_y2 = anticommutator(obj.x2, variate);
    
    % Second term of ∇({x_1^4, x_2^2}) -> 0.5 * {x_1^4, {x_2, y_2}}
    expr = 0.5 * anticommutator(obj.x1_pow4, ac_x2_y2);
    
    % First term of ∇({x_1^2, x_2^4}) -> 0.5 * {x_1^2, {x_2^2, {x_2, y_2}}}
    expr = expr + ...
        0.5 * anticommutator(obj.x1_pow2, ...
                             anticommutator(obj.x2_pow2, ac_x2_y2));
    
    % Second term of ∇(x_1^2 x_2^2 x_3^2) -> -1.5 * x_1^2 * {x_2, y_2} x_3^2 
    expr = expr - 1.5 * obj.x1_pow2 * ac_x2_y2 * obj.x3_pow2;
    
    % Second term of ∇(x_3^2 x_2^2 x_1^2) -> -1.5 * x_3^2 * {x_2, y_2} x_1^2 
    expr = expr - 1.5 * obj.x3_pow2 * ac_x2_y2 * obj.x1_pow2;
    
end

function expr = nabla_objective_y3(obj, variate) 
%NABLA_OBJECTIVE_Y3 y_3 specialization of derivative.

    % {x_3, y_3}:
    ac_x3_y3 = anticommutator(obj.x3, variate);
    
    % First term of ∇(x_3^6) -> {x_3^4, {x_3, y_3}}
    expr = anticommutator(obj.x3_pow4, ac_x3_y3);
    
    % Second term of ∇(x_3^6) -> x_3^2 {x_3, y_3} x_3^2
    expr = expr + obj.x3_pow2 * ac_x3_y3 * obj.x3_pow2;
    
    % First term of ∇(x_1^2 x_2^2 x_3^2) -> -1.5 * x_1^2 x_2^2 {x_3, y_3}
    expr = expr - 1.5 * obj.x1_pow2 * obj.x2_pow2 * ac_x3_y3;
    
    % Third term of ∇(x_3^2 x_2^2 x_1^2) -> -1.5 * {x_3, y_3} x_2^2 x_1^2 
    expr = expr - 1.5 * ac_x3_y3 * obj.x2_pow2 * obj.x1_pow2;

end