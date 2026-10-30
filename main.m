% MAIN  B1 Numerical Algorithms mini-project: Parts A, B and C.

clear; clc; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Setup: parameters and noisy "sensor" signal
rng(1);                                    % required for reproducibility; fixes randn and randperm
P = params();
y = P.ytrue + P.sigma * randn(P.N, 1);     % first use of the global stream (before Part B's split)

%% Part A - Numerical Differentiation
A = partA_differentiation(P, y);

%% Part B - Polynomial Approximation
B = partB_poly(P, y);

%% Part C - Numerical ODE Solutions
C = partC_ODE(P);

%% Output
print_results(P, A, B, C);
plot_figures(P, A, B, C);
