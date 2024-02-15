%% compare.m
% Produces the values featured in arXiv:2311.18707.
%

%% Parameters
d = 3;
delta = 0.1;
objective_tensor = [0, 1,  1, 0; ...
                    1, 1,  1, 1; ...
                    1, 1,  1, -1; ...
                    0, 1, -1,  0];
essential_epsilon = 1e-3;
mm_level = 2;
lm_level = 2;
kkt_level = 3; % Ideally: mm_level * 2 - 1
verbose = 1;   % Set to 0 (quiet), 1 (some output) or 2 (lots of output).

%% Solve
problem = MFCQExample(d, delta, objective_tensor, verbose);

% No operator ncKKT, no state-optimality
without_kkt = problem.solve_without_kkt(mm_level, lm_level, false);

% No operator ncKKT, with state-optimality
without_kkt_so = problem.solve_without_kkt(mm_level, lm_level, true);

% Essential ncKKT, no state-optimality
essential_kkt = problem.solve_essential_kkt(mm_level, lm_level, ...
                                            kkt_level, ...
                                            essential_epsilon, false);
                                        
% Essential ncKKT, with state-optimality
essential_kkt_so = problem.solve_essential_kkt(mm_level, lm_level, ...
                                               kkt_level, ...
                                               essential_epsilon, true);

% Weak ncKKT, no state-optimality
weak_kkt = problem.solve_weak_kkt(mm_level, lm_level, kkt_level, false);

% Weak ncKKT, with state-optimality
weak_kkt_so = problem.solve_weak_kkt(mm_level, lm_level, kkt_level, true);

% Normed ncKKT, no state-optimality
normed_kkt = problem.solve_normed_kkt(mm_level, lm_level, kkt_level, false);

% Normed ncKKT, with state-optimality
normed_kkt_so = problem.solve_normed_kkt(mm_level, lm_level, kkt_level, true);

%% Display summary of results
fprintf("\nResults for delta = %g, MM = %d, LM = %d, KKT = %d:\n", ...
        delta, mm_level, lm_level, kkt_level);
fprintf("Without ncKKT:\n%.10g\n\n", without_kkt);
fprintf("Without ncKKT, with state-optimality:\n%.10g\n\n", without_kkt_so);
fprintf("Essential ncKKT (epsilon = %g):\n%.10g\n\n", ...
        essential_epsilon, essential_kkt);
fprintf("Essential ncKKT (epsilon = %g), with state-optimality:\n%.10g\n\n", ...
        essential_epsilon, essential_kkt_so);
fprintf("Weak ncKKT:\n%.10g\n\n", weak_kkt);
fprintf("Weak ncKKT, with state-optimality:\n%.10g\n\n", weak_kkt_so);
fprintf("Normed ncKKT:\n%.10g\n\n", normed_kkt);
fprintf("Normed ncKKT, with state-optimality:\n%.10g\n\n", normed_kkt_so);