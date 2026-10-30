function [t, Y] = rk4_solve(f, Y0, T, N)
% RK4_SOLVE  4th-order Runge-Kutta for Y' = f(t, Y) using N steps.
%   update is:
%   Y(n+1) = Y(n) + h/6 * (k1 + 2*k2 + 2*k3 + k4).
%
%   Inputs:
%     f  - [function handle] f(t, Y) -> [d x 1] derivative; t scalar, Y [d x 1]
%     Y0 - [d x 1] initial state
%     T  - [scalar] final time
%     N  - [scalar int] number of steps
%   Outputs:
%     t  - [1 x (N+1)] time points
%     Y  - [d x (N+1)] state at each time point

    h = T / N;
    t = (0:N) * h;
    
    Y = zeros(numel(Y0), N + 1);
    Y(:, 1) = Y0(:);
    
    for n = 1:N
        k1 = f(t(n),       Y(:, n));
        k2 = f(t(n) + h/2, Y(:, n) + h/2 * k1);
        k3 = f(t(n) + h/2, Y(:, n) + h/2 * k2);
        k4 = f(t(n) + h,   Y(:, n) + h   * k3);
    
        Y(:, n+1) = Y(:, n) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end
