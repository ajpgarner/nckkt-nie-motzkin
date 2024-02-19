%% example_534.m
%

%% Parameters
delta = 0.01;
min_radius = 1.0;
max_radius = 10.0;
espilon = 1e-4; % For essential-ncKKT
verbosity = 1;

%% Set-up
solver = NCNieMotzkin(delta, min_radius, max_radius, verbosity);

%% Solve:
nckkt_534 = solver.solve_kkt(5, 3, 4, espilon, false);
fprintf("KKT,\tMM = %d,\tLM = %d,\tKKT = %d:\t%.8g\n", 5, 3, 4, nckkt_534);
fprintf("\n");

