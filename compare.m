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

%% NPA Solves:
npa_33 = solver.solve_without_kkt(3, 3, false);
npa_43 = solver.solve_without_kkt(4, 3, false);

%% ncKKT Solves:
nckkt_331 = solver.solve_kkt(3, 3, 1, espilon, false);
nckkt_431 = solver.solve_kkt(4, 3, 1, espilon, false);
nckkt_432 = solver.solve_kkt(4, 3, 2, espilon, false);


%% Display summary of results
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 3, 3, npa_33);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 4, 3, npa_43);
fprintf("\n");

fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d,\t%epsilon=%g:\t%.8g\n", 3, 3, 1, nckkt_331);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d,\t%epsilon=%g:\t%.8g\n", 4, 3, 1, nckkt_431);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d,\t%epsilon=%g:\t%.8g\n", 4, 3, 2, nckkt_432);

fprintf("\n");

