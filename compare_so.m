%% compare.m
%

%% Parameters
delta = 0.01;
min_radius = 1.0;
max_radius = 7.5;
epsilon = 1e-4; % For essential-ncKKT
verbosity = 1;

%% Set-up
solver = NCNieMotzkin(delta, min_radius, max_radius, verbosity);

%% Solve:
so_npa_33 = solver.solve_without_kkt(3, 3, 4);
so_npa_43 = solver.solve_without_kkt(4, 3, 4);
so_nckkt_331 = solver.solve_kkt(3, 3, 1, epsilon, 4);
so_nckkt_433 = solver.solve_kkt(4, 3, 3, epsilon, 4);
so_npa_53 = solver.solve_without_kkt(5, 3, 4);


%% Display summary of results
fprintf("No KKT,\tMM = %d,\tLM = %d,\tSO = %d:\t%.8g\n", 3, 3, 4, so_npa_33);
fprintf("No KKT,\tMM = %d,\tLM = %d,\tSO = %d:\t%.8g\n", 4, 3, 4, so_npa_43);
fprintf("No KKT,\tMM = %d,\tLM = %d,\tSO = %d:\t%.8g\n", 5, 3, 4, so_npa_53);
fprintf("\n");

fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d,\tSO = %d,\tepsilon=%g:\t%.8g\n", 3, 3, 1, 4, epsilon, so_nckkt_331);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d,\tSO = %d,\tepsilon=%g:\t%.8g\n", 4, 3, 3, 4, epsilon, so_nckkt_433);

fprintf("\n");

