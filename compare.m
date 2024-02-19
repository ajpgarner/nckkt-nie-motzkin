%% compare.m
%

%% Parameters
delta = 0.01;
min_radius = 1.0;
max_radius = 7.5;
espilon = 1e-4; % For essential-ncKKT
verbosity = 1;

%% Set-up
solver = NCNieMotzkin(delta, min_radius, max_radius, verbosity);

%% Solve:
npa_33 = solver.solve_without_kkt(3, 3, false);
npa_43 = solver.solve_without_kkt(4, 3, false);
nckkt_331 = solver.solve_kkt(3, 3, 1, espilon, false);
nckkt_433 = solver.solve_kkt(4, 3, 3, espilon, false);
npa_53 = solver.solve_without_kkt(5, 3, false);


%% Display summary of results
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 3, 3, npa_33);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 4, 3, npa_43);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 5, 3, npa_53);
fprintf("\n");

fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d,\t%epsilon=%g:\t%.8g\n", 3, 3, 1, epsilon, nckkt_331);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d,\t%epsilon=%g:\t%.8g\n", 4, 3, 3, epsilon, nckkt_433);

fprintf("\n");

