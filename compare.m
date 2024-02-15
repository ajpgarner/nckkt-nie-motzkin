%% compare.m
%

%% Parameters
delta = 0.01;
min_radius = 1.0;
max_radius = 10.0;
espilon = 1e-4; % For essential-ncKKT
verbosity = 1;

%% Set-up
solver = NCNie(delta, min_radius, max_radius, verbosity);

%% NPA Solves:
npa_32 = solver.solve_without_kkt(3, 2, false);
npa_33 = solver.solve_without_kkt(3, 3, false);
npa_43 = solver.solve_without_kkt(4, 3, false);
npa_44 = solver.solve_without_kkt(4, 4, false);
npa_54 = solver.solve_without_kkt(5, 4, false);

%% ncKKT Solves:
nckkt_321 = solver.solve_kkt(3, 2, 1, espilon, false);
nckkt_332 = solver.solve_kkt(3, 3, 2, espilon, false);
nckkt_433 = solver.solve_kkt(4, 3, 3, espilon, false);
nckkt_443 = solver.solve_kkt(4, 4, 3, espilon, false);
nckkt_544 = solver.solve_kkt(5, 4, 4, espilon, false);
nckkt_545 = solver.solve_kkt(5, 4, 5, espilon, false);


%% Display summary of results
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 3, 2, npa_32);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 3, 3, npa_3);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 4, 3, npa_43);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 4, 4, npa_44);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 5, 4, npa_54);
fprintf("\n");

fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 3, 2, 1, nckkt_321);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 3, 3, 2, nckkt_332);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 4, 3, 3, nckkt_433);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 4, 4, 3, nckkt_443);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 5, 4, 4, nckkt_544);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 5, 4, 5, nckkt_545);
fprintf("\n");

