%% compare.m
%

%% Parameters
delta = 0.01;
min_radius = 1.0;
max_radius = 10.0;
espilon = 1e-3: % For essential-ncKKT.
verbosity = 2;

%% Set-up
solver = NCNie(delta, min_radius, max_radius, verbosity);

%% NPA Solves:
npa_32 = solver.solve_without_kkt(3, 2, false);
npa_33 = solver.solve_without_kkt(3, 3, false);
npa_43 = solver.solve_without_kkt(4, 3, false);
npa_44 = solver.solve_without_kkt(4, 4, false);
npa_54 = solver.solve_without_kkt(5, 4, false);

%% ncKKT Solves:
nckkt_325 = solver.solve_without_kkt(3, 2, 5, espilon, false);
nckkt_335 = solver.solve_without_kkt(3, 3, 5, espilon, false);
nckkt_437 = solver.solve_without_kkt(4, 3, 7, espilon, false);
nckkt_447 = solver.solve_without_kkt(4, 4, 7, espilon, false);


%% Display summary of results
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 3, 2, npa_32);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 3, 3, npa_3);
fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 4, 3, npa_43);
%fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 4, 4, npa_44);
%fprintf("No KKT,\tMM = %d,\tLM = %d:\t%.8g\n", 5, 4, npa_54);
fprintf("\n");

fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 3, 2, 5, nckkt_325);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 3, 3, 5, nckkt_335);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 4, 3, 7, nckkt_437);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 4, 4, 7, nckkt_447);
fprintf("\n");

